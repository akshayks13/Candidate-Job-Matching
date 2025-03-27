#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <cctype>
#include <limits>
#include <cmath>
#include <nlohmann/json.hpp>
#include <iomanip>

const double MATCH_THRESHOLD = 100.0;

using json = nlohmann::json;
using namespace std;

// Certificate structure
struct Certificate {
    string title;
    string issuer;
};

// Candidate structure
struct Candidate {
    int id;
    string name;
    string resume;           
    vector<string> skills;   
    double cgpa;
    int experience;
    vector<Certificate> certificates;
    int expectedSalary;
    string preferredLocation;
    vector<string> softSkills;
};

// Job structure
struct Job {
    int id;
    string title;
    string description;            
    vector<string> requiredSkills; 
    double min_cgpa;
    int min_experience;
    vector<string> requiredCertifications;
    int offeredSalary;
    string location;
    vector<string> requiredSoftSkills;
};

// MatchResult structure
struct MatchResult {
    Candidate candidate;
    Job job;
    double cost;
    bool suitable; // true if cost <= MATCH_THRESHOLD
};

// Read JSON from file
json readJSON(const string &filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        exit(1);
    }
    json j;
    file >> j;
    return j;
}

// Convert a string to lowercase
string toLower(const string &s) {
    string result = s;
    transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

// Tokenize a string into a set of words
set<string> tokenize(const string &s) {
    // Define a set of common stopwords to ignore.
    const set<string> STOPWORDS = {
        "the", "and", "a", "an", "of", "in", "on", "for", "with", "to", "at", "by", "from", "that", "this", "is", "it", "as"
    };
    
    set<string> tokens;
    string word;
    for (char c : s) {
        if (isalnum(c)) {
            word.push_back(tolower(c));
        } else if (!word.empty()) {
            // Only insert if word is not a stopword.
            if (STOPWORDS.find(word) == STOPWORDS.end())
                tokens.insert(word);
            word.clear();
        }
    }
    // Check any trailing word.
    if (!word.empty() && STOPWORDS.find(word) == STOPWORDS.end())
        tokens.insert(word);
    return tokens;
}

// Jaccard similarity calculation between two token sets
double jaccardSimilarity(const set<string> &a, const set<string> &b) {
    int intersection = 0;
    for (const auto &token : a)
        if (b.find(token) != b.end())
            intersection++;
    int unionSize = a.size() + b.size() - intersection;
    return (unionSize == 0) ? 0.0 : static_cast<double>(intersection) / unionSize;
}

set<string> removeCommonTokens(const set<string>& tokens) {
    // Define a set of common tokens to ignore.
    static const set<string> commonTokens = {"certified", "certificate", "certification", "professional"};
    set<string> filtered;
    for (const string &token : tokens) {
        if (commonTokens.find(token) == commonTokens.end()) {
            filtered.insert(token);
        }
    }
    return filtered;
}

// Check if a certificate is relevant to the job by comparing its title with required certification strings.
// We use Jaccard similarity on token sets and consider it relevant if the similarity exceeds a threshold.
bool isCertificateRelevant(const Certificate &cert, const Job &j, double threshold = 0.5) {
    // Tokenize the candidate certificate title.
    set<string> certTokens = tokenize(cert.title);
    // For each required certification from the job...
    for (const string &reqCert : j.requiredCertifications) {
        set<string> reqTokens = tokenize(reqCert);
        reqTokens = removeCommonTokens(reqTokens);
        double sim = jaccardSimilarity(certTokens, reqTokens);
        if (sim >= threshold) {
            return true;
        }
    }
    return false;
}

// Calculate the certification cost based on the relevance of candidate's certificates to the job requirements.
int calculateCertificationCost(const Candidate &c, const Job &j) {
    // If the job doesn't require certifications, return zero cost.
    if (j.requiredCertifications.empty()) {
        return 0;
    }
    
    int relevantCount = 0;
    // Check each candidate certificate.
    for (const Certificate &cert : c.certificates) {
        if (isCertificateRelevant(cert, j)) {
            relevantCount++;
        }
    }
    
    int requiredCount = j.requiredCertifications.size();
    int cost = 0;
    if (relevantCount < requiredCount) {
        // Penalize for each missing relevant certificate.
        cost += (requiredCount - relevantCount) * 5;
    } else if (relevantCount > requiredCount) {
        // Provide a bonus for extra relevant certificates.
        cost -= (relevantCount - requiredCount) * 3;
    }
    return cost;
}

// Check if one string is a substring of the other.
bool isSubstringMatch(const string &s1, const string &s2) {
    string s1Lower = toLower(s1);
    string s2Lower = toLower(s2);
    return (s1Lower.find(s2Lower) != string::npos) || (s2Lower.find(s1Lower) != string::npos);
}

// Check if two skill strings match using a combination of substring matching and Jaccard similarity.
bool isSkillMatch(const string &candidateSkill, const string &jobSkill, double jaccardThreshold = 0.5) {
    // First, try substring matching.
    if (isSubstringMatch(candidateSkill, jobSkill))
        return true;
    
    // Otherwise, compute Jaccard similarity on token sets.
    set<string> candTokens = tokenize(candidateSkill);
    set<string> jobTokens = tokenize(jobSkill);
    double sim = jaccardSimilarity(candTokens, jobTokens);
    return sim >= jaccardThreshold;
}

// ----- Structured Skill Matching using Jaccard Similarity -----
// This function compares candidate.skills vector with job.requiredSkills.
double calculateSkillCost(const Candidate &c, const Job &j) {
    int matchingRequiredSkillCount = 0;
    
    // For each required skill, see if any candidate skill matches using hybrid matching.
    for (const string &reqSkill : j.requiredSkills) {
        bool foundMatch = false;
        for (const string &candSkill : c.skills) {
            if (isSkillMatch(candSkill, reqSkill)) {
                foundMatch = true;
                break;  // Found a match for this required skill.
            }
        }
        if (foundMatch)
            matchingRequiredSkillCount++;
    }
    
    // Compute penalty based on missing required skills.
    int missing = j.requiredSkills.size() - matchingRequiredSkillCount;
    double penalty = missing * 10; 
    return penalty;
}

// ----- Soft Skill Matching using Jaccard Similarity -----
// The function compares candidates soft skills with the jobs requirements
double calculateSoftSkillCost(const Candidate &c, const Job &j) {
    if (j.requiredSoftSkills.empty()) return 0;
    
    set<string> candidateSoftSet, jobSoftSet;
    // Use candidate.softSkills if available; otherwise, extract from resume.
    vector<string> candidateSoft = c.softSkills.empty() ? vector<string>() : c.softSkills;

    for (const string &s : candidateSoft)
        candidateSoftSet.insert(toLower(s));
    for (const string &s : j.requiredSoftSkills)
        jobSoftSet.insert(toLower(s));

    if (candidateSoftSet.empty()) {
        // For example, a fixed penalty of 10 points if no soft skills are provided.
        return 10;
    }
    double sim = jaccardSimilarity(candidateSoftSet, jobSoftSet);
    double scale = 10; 
    return -(sim * scale); 
}

// Structured cost function based on skills, experience, salary, and location
int calculateStructuredCost(const Candidate &c, const Job &j) {
    int cost = 0;

    int matchingRequiredSkillCount = 0;
    for (const string &reqSkill : j.requiredSkills) {
        if (find(c.skills.begin(), c.skills.end(), reqSkill) != c.skills.end()) {
            matchingRequiredSkillCount++;
        }
    }
    
    // Enforce a minimum threshold: if the candidate doesn't have any required skill, return a very high cost.
    if (matchingRequiredSkillCount == 0) {
        return 1e6; 
    }
    
    // Skill mismatch penalty: start with a penalty equal to required skills, reduce for each matching skill.
    int skillPenalty = j.requiredSkills.size();
    for (const string &skill : c.skills) {
        if (find(j.requiredSkills.begin(), j.requiredSkills.end(), skill) != j.requiredSkills.end()) {
            skillPenalty--; // one less missing skill
        }
    }
    cost += skillPenalty * 10;

    // Experience penalty: penalize if candidate has less than required experience.
    if (c.experience < j.min_experience) {
        cost += (j.min_experience - c.experience) * 5;
    }

    // Salary penalty: if candidate's expected salary is higher than offered.
    if (c.expectedSalary > j.offeredSalary) {
        cost += (c.expectedSalary - j.offeredSalary) / 1000; 
    }

    // Location penalty: fixed penalty for location mismatch.
    if (j.location != "Remote" && c.preferredLocation != j.location) {
        cost += 15;
    }
    
    // Certification penalty: call the improved certification matching function.
    cost += calculateCertificationCost(c, j);
    cost += calculateSoftSkillCost(c, j);

    return cost;
}

// Compute Smith Waterman similarity between two strings.
// Returns a normalized similarity value between 0 and 1.
double smithWatermanSimilarity(const string &s1, const string &s2) {
    const int matchScore = 2;
    const int mismatchPenalty = -1;
    const int gapPenalty = -1;
    
    int m = s1.size();
    int n = s2.size();
    
    // Create a (m+1) x (n+1) scoring matrix initialized to 0.
    vector<vector<int>> score(m + 1, vector<int>(n + 1, 0));
    int maxScore = 0;
    
    // Fill in the scoring matrix.
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int scoreDiag = score[i-1][j-1] + (tolower(s1[i-1]) == tolower(s2[j-1]) ? matchScore : mismatchPenalty);
            int scoreUp = score[i-1][j] + gapPenalty;
            int scoreLeft = score[i][j-1] + gapPenalty;
            score[i][j] = max(0, max(scoreDiag, max(scoreUp, scoreLeft)));
            if (score[i][j] > maxScore) {
                maxScore = score[i][j];
            }
        }
    }
    
    // Normalize the score. The maximum possible score is (min(m, n) * matchScore).
    int maxPossible = min(m, n) * matchScore;
    return (maxPossible > 0) ? static_cast<double>(maxScore) / maxPossible : 0.0;
}

// Hybrid cost function combining text-based and structured matching.
// Uses Smith Waterman similarity for Resume-Job Description matching.
double hybridCost(const Candidate &c, const Job &j) {
    // Compute text-based similarity using Smith Waterman.
    double textSim = smithWatermanSimilarity(c.resume, j.description);

    double costText = (1.0 - textSim) * 50;  

    // Compute structured cost.
    int costStructured = calculateStructuredCost(c, j);

    // Combine the two costs using weights.
    double alpha = 0.6;  // Weight for text-based(Resume-Job Description) cost
    double beta  = 0.4;  // Weight for structured(Other fields) cost

    return alpha * costText + beta * costStructured;
}


// Hungarian Algorithm for optimal assignment (minimization)
vector<int> hungarian(const vector<vector<double>> &cost) {
    int n = cost.size(), m = cost[0].size();
    int size = max(n, m);
    vector<vector<double>> a(size, vector<double>(size, 0));
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            a[i][j] = (i < n && j < m) ? cost[i][j] : 0;

    vector<double> u(size+1), v(size+1);
    vector<int> p(size+1), way(size+1);

    for (int i = 1; i <= size; i++) {
        p[0] = i;
        double inf = numeric_limits<double>::infinity();
        vector<double> minv(size+1, inf);
        vector<bool> used(size+1, false);
        int j0 = 0;
        do {
            used[j0] = true;
            int i0 = p[j0], j1 = 0;
            double delta = inf;
            for (int j = 1; j <= size; j++) {
                if (!used[j]) {
                    double cur = a[i0-1][j-1] - u[i0] - v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            for (int j = 0; j <= size; j++) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }

    vector<int> assignment(n, -1);
    for (int j = 1; j <= size; j++)
        if (p[j] <= n && j <= m) assignment[p[j]-1] = j-1;
    return assignment;
}

// Function to display the match analysis between a candidate and a job.
// Will say about why the candidate is suitable or not for the job.
void displayMatchAnalysis(const Candidate &cand, const Job &job, double cost, bool suitable) {
    cout << "\nMatch Analysis for Candidate: " << cand.name << " (ID: " << cand.id << ") and Job: " << job.title << " (ID: " << job.id << ")\n";
    cout << string(80, '=') << endl;
    
    // Display table headers
    cout << left << setw(25) << "Criteria" 
         << setw(25) << "Candidate" 
         << setw(25) << "Job Requirement" 
         << setw(15) << "Match" << endl;
    cout << string(80, '-') << endl;
    
    // Skills Matching
    int matchCount = 0;
    for (const string &reqSkill : job.requiredSkills) {
        for (const string &candSkill : cand.skills) {
            if (isSkillMatch(candSkill, reqSkill)) {
                matchCount++;
                break;
            }
        }
    }
    // 'Yes' only if ALL required skills are matched.
    bool hasSkillMatch = (matchCount == job.requiredSkills.size());
    cout << left << setw(25) << "Skills Matched"
        << setw(25) << to_string(matchCount) + " / " + to_string(job.requiredSkills.size())
        << setw(25) << "Required: " + to_string(job.requiredSkills.size()) 
        << setw(15) << (hasSkillMatch ? "Yes" : "No") << endl;

    // Experience Matching
    bool expMatch = cand.experience >= job.min_experience;
    cout << left << setw(25) << "Experience"
         << setw(25) << to_string(cand.experience) + " years"
         << setw(25) << "Min: " + to_string(job.min_experience) + " years"
         << setw(15) << (expMatch ? "Yes" : "No") << endl;
    
    // Salary Expectation Check
    bool salaryMatch = cand.expectedSalary <= job.offeredSalary;
    cout << left << setw(25) << "Salary Expectation"
         << setw(25) << to_string(cand.expectedSalary)
         << setw(25) << "Max: " + to_string(job.offeredSalary)
         << setw(15) << (salaryMatch ? "Yes" : "No") << endl;
    
    // Location Match
    bool locationMatch = (job.location == "Remote" || cand.preferredLocation == job.location);
    cout << left << setw(25) << "Location Match"
         << setw(25) << cand.preferredLocation
         << setw(25) << job.location
         << setw(15) << (locationMatch ? "Yes" : "No") << endl;
    
    // Certifications Match
    int certMatchCount = 0;
    for (const Certificate &cert : cand.certificates) {
        if (isCertificateRelevant(cert, job))
            certMatchCount++;
    }
    // 'Yes' only if the candidate meets or exceeds all required certifications.
    bool certMatch = (certMatchCount >= job.requiredCertifications.size());
    cout << left << setw(25) << "Certifications Matched"
        << setw(25) << to_string(certMatchCount)
        << setw(25) << "Required: " + to_string(job.requiredCertifications.size())
        << setw(15) << (certMatch ? "Yes" : "No") << endl;
    
    // Soft Skills Similarity
    set<string> candidateSoftSet, jobSoftSet;
    for (const string &s : cand.softSkills) candidateSoftSet.insert(toLower(s));
    for (const string &s : job.requiredSoftSkills) jobSoftSet.insert(toLower(s));
    double softSim = jaccardSimilarity(candidateSoftSet, jobSoftSet);
    cout << left << setw(25) << "Soft Skills Similarity"
         << setw(25) << to_string(softSim)
         << setw(25) << "Threshold: 0.5"
         << setw(15) << (softSim >= 0.5 ? "Yes" : "No") << endl;
    
    // Resume-Job Description Similarity
    double textSim = smithWatermanSimilarity(cand.resume, job.description);
    cout << left << setw(25) << "Resume-Job Similarity"
         << setw(25) << to_string(textSim)
         << setw(25) << "Threshold: 0.5"
         << setw(15) << (textSim >= 0.5 ? "Yes" : "No") << endl;
    
    cout << string(80, '=') << endl;
    cout << left << setw(25) << "Hybrid Cost" 
         << setw(25) << to_string(cost) 
         << setw(25) << "Threshold: " + to_string(MATCH_THRESHOLD) 
         << setw(15) << (suitable ? "Suitable" : "Not Suitable") << endl;
    
    cout << string(80, '=') << endl;
}

// Function to generate a match between a candidate and a job.
MatchResult generateMatch(const Candidate &cand, const Job &job) {
    double cost = hybridCost(cand, job);
    bool suitable = (cost <= MATCH_THRESHOLD);
    return {cand, job, cost, suitable};
}

// Function to search for candidates that have all of the given skills.
void searchCandidatesByMultipleSkills(const vector<Candidate>& candidates) {
    int numSkills;
    cout << "How many skills do you want to enter? ";
    cin >> numSkills;
    cin.ignore(); 

    vector<string> requiredSkills;
    for (int i = 0; i < numSkills; i++) {
        string skill;
        cout << "Enter skill #" << (i + 1) << ": ";
        getline(cin, skill);
        requiredSkills.push_back(toLower(skill));
    }
    
    cout << "\nCandidates that have all these skills:" << endl;
    bool foundCandidate = false;
    for (const auto &cand : candidates) {
        vector<string> candSkills;
        for (const string &s : cand.skills) {
            candSkills.push_back(toLower(s));
        }
        
        bool hasAll = true;
        for (const auto &reqSkill : requiredSkills) {
            bool skillFound = false;
            for (const auto &candSkill : candSkills) {
                if (candSkill == reqSkill) {
                    skillFound = true;
                    break;
                }
            }
            if (!skillFound) {
                hasAll = false;
                break;
            }
        }
        
        if (hasAll) {
            cout << "Candidate ID: " << cand.id << ", Name: " << cand.name << endl;
            foundCandidate = true;
        }
    }
    
    if (!foundCandidate) {
        cout << "No candidate found with all the entered skills." << endl;
    }
}

void displayMenu() {
    cout << "\nJob Matching System" << endl;
    cout << "1. Candidate List" << endl;
    cout << "2. Job List" << endl;
    cout << "3. List most suitable jobs for each candidate" << endl;
    cout << "4. List all suitable jobs for a given Candidate (sorted by suitability)" << endl;
    cout << "5. List all suitable Candidates for a given job (sorted by suitability)" << endl;
    cout << "6. Match Report" << endl;
    cout << "7. Search by skill" << endl;
    cout << "8. Exit" << endl;
    cout << "Enter your choice: ";
}

int main() {
    
    // Read JSON input
    json data = readJSON("input.json");

    vector<Candidate> candidates;
    vector<Job> jobs;

    // Parse candidate entries
    for (auto &c : data["candidates"]) {
        Candidate cand;
        cand.id = c["id"];
        cand.name = c["name"];
        cand.resume = c["resume"];
        for (auto &skill : c["skills"])
            cand.skills.push_back(skill);
        cand.cgpa = c["cgpa"];
        cand.experience = c["experience"];
        for (auto &cert : c["certificates"]) {
            Certificate certificate;
            certificate.title = cert["title"];
            certificate.issuer = cert["issuer"];
            cand.certificates.push_back(certificate);
        }
        cand.expectedSalary = c["expectedSalary"];
        cand.preferredLocation = c["preferredLocation"];
        for(auto &soft : c["softSkills"]) {
            cand.softSkills.push_back(soft);
        }
        candidates.push_back(cand);
    }

    // Parse job entries
    for (auto &j : data["jobs"]) {
        Job job;
        job.id = j["id"];
        job.title = j["title"];
        job.description = j["description"];
        for (auto &req : j["requiredSkills"])
            job.requiredSkills.push_back(req);
        job.min_cgpa = j["min_cgpa"];
        job.min_experience = j["min_experience"];
        for (auto &reqCert : j["requiredCertifications"])
            job.requiredCertifications.push_back(reqCert);
        job.offeredSalary = j["offeredSalary"];
        job.location = j["location"];
        for(auto &soft : j["requiredSoftSkills"]) {
            job.requiredSoftSkills.push_back(soft);
        }
        jobs.push_back(job);
    }

    // Build cost matrix using the hybrid cost function.
    int numCandidates = candidates.size();
    int numJobs = jobs.size();
    int size = max(numCandidates, numJobs);
    vector<vector<double>> costMatrix(numCandidates, vector<double>(numJobs, 0));
    for (int i = 0; i < numCandidates; i++) {
        for (int j = 0; j < numJobs; j++) {
            costMatrix[i][j] = hybridCost(candidates[i], jobs[j]);
        }
    }

    vector<vector<double>> paddedCost(size, vector<double>(size, 0));
    // Copy the existing cost values.
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i < numCandidates && j < numJobs) {
                paddedCost[i][j] = costMatrix[i][j];
            } else {
                // Assign a high cost to dummy entries.
                paddedCost[i][j] = 1e6; 
            }
        }
    }

    // Run the Hungarian algorithm to get optimal assignment.
    vector<int> assignment = hungarian(paddedCost);

    cout << "Cost Matrix:" << endl;
    for (const auto &row : paddedCost) {
        for (double val : row)
            cout << val << "\t";
        cout << endl;
    }
    cout << endl;

    int choice;
    while (true) {
        displayMenu();
        cin >> choice;

        if (choice == 1) {
            // Option 1: List all Candidate names and their IDs.
            cout << "\nCandidate List:" << endl;
            for (const auto &cand : candidates) {
                cout << "ID: " << cand.id << " - Name: " << cand.name << endl;
            }
        }else if (choice == 2) {
            // Option 2: List all job names and their IDs.
            cout << "\nJob List:" << endl;
            for (const auto &job : jobs) {
                cout << "ID: " << job.id << " - Title: " << job.title << endl;
            }
        }else if (choice == 3) {
            // Print assignment results in tabular form.
            cout << "Optimal Assignment:" << endl;
            
            // Print header row.
            cout << left << setw(12) << "Cand. ID" 
                 << setw(25) << "Candidate Name" 
                 << setw(10) << "Job ID" 
                 << setw(30) << "Job Title" 
                 << setw(15) << "Hybrid Cost" << endl;
            
            cout << string(92, '-') << endl;
        
            // Iterate over each candidate.
            for (int i = 0; i < numCandidates; i++) {
                if (assignment[i] < numJobs && paddedCost[i][assignment[i]] <= MATCH_THRESHOLD) {
                    cout << left << setw(12) << candidates[i].id 
                         << setw(25) << candidates[i].name 
                         << setw(10) << jobs[assignment[i]].id 
                         << setw(30) << jobs[assignment[i]].title 
                         << setw(15) << paddedCost[i][assignment[i]] << endl;
                } else {
                    cout << left << setw(12) << candidates[i].id 
                         << setw(25) << candidates[i].name 
                         << setw(10) << "N/A" 
                         << setw(30) << "No job selected" 
                         << setw(15) << "N/A" << endl;
                }
            }
        }else if (choice == 4) {
            // Option: Ask for a Candidate's name and ID, then list jobs categorized by suitability.
            string empName;
            int empId;
            cout << "Enter Candidate name: ";
            cin.ignore();
            getline(cin, empName);
            cout << "Enter Candidate ID: ";
            cin >> empId;
            
            bool found = false;
            Candidate targetCandidate;
            for (const auto &cand : candidates) {
                if (cand.name == empName && cand.id == empId) {
                    targetCandidate = cand;
                    found = true;
                    break;
                }
            }
            
            if (!found) {
                cout << "Candidate not found!" << endl;
            } else {
                vector<pair<double, Job>> suitableRanking;
                vector<pair<double, Job>> notSuitableRanking;
                
                // For each job, calculate the hybrid cost for this candidate.
                for (const auto &job : jobs) {
                    double cost = hybridCost(targetCandidate, job);
                    if (cost <= MATCH_THRESHOLD)
                        suitableRanking.push_back({cost, job});
                    else
                        notSuitableRanking.push_back({cost, job});
                }
                
                // Sort both categories by ascending cost.
                sort(suitableRanking.begin(), suitableRanking.end(), [](const pair<double, Job>& a, const pair<double, Job>& b) {
                    return a.first < b.first;
                });
                sort(notSuitableRanking.begin(), notSuitableRanking.end(), [](const pair<double, Job>& a, const pair<double, Job>& b) {
                    return a.first < b.first;
                });
                
                // Display suitable jobs.
                if (suitableRanking.empty()) {
                    cout << "\nNo job is considered suitable for " << targetCandidate.name << " (Hybrid Cost <= " << MATCH_THRESHOLD << ")." << endl;
                } else {
                    cout << "\nJobs suitable for " << targetCandidate.name << " (from most to least suitable):" << endl;
                    for (const auto &entry : suitableRanking) {
                        cout << "Job ID: " << entry.second.id 
                             << " - Job: " << entry.second.title 
                             << " (Hybrid Cost: " << entry.first << ")" << endl;
                    }
                }
                
                // Display the not suitable jobs.
                if (!notSuitableRanking.empty()) {
                    cout << "\nJobs not suitable for " << targetCandidate.name << " (Hybrid Cost > " << MATCH_THRESHOLD << "):" << endl;
                    for (const auto &entry : notSuitableRanking) {
                        cout << "Job ID: " << entry.second.id 
                             << " - Job: " << entry.second.title 
                             << " (Hybrid Cost: " << entry.first << ")" << endl;
                    }
                }
            }
        }else if (choice == 5) {
            // Option 5: Find the most suitable Candidates for a given job.
            string jobTitle;
            int jobId;
            
            // Prompt for both job title and job ID.
            cout << "Enter job title: ";
            cin.ignore(); 
            getline(cin, jobTitle);
            cout << "Enter job ID: ";
            cin >> jobId;
            
            bool found = false;
            Job targetJob;
            // Search for the job using both job title and job ID.
            for (const auto &job : jobs) {
                if (toLower(job.title) == toLower(jobTitle) && job.id == jobId) {
                    targetJob = job;
                    found = true;
                    break;
                }
            }
            
            if (!found) {
                cout << "Job not found!" << endl;
            } else {
                vector<pair<double, Candidate>> suitableCandidates;
                vector<pair<double, Candidate>> notSuitableCandidates;
                
                // For each candidate, calculate the hybrid cost for this job.
                for (const auto &cand : candidates) {
                    double cost = hybridCost(cand, targetJob);
                    if (cost <= MATCH_THRESHOLD)
                        suitableCandidates.push_back({cost, cand});
                    else
                        notSuitableCandidates.push_back({cost, cand});
                }
                
                // Sort both lists by ascending cost (lower cost means more suitable).
                sort(suitableCandidates.begin(), suitableCandidates.end(), [](const pair<double, Candidate>& a, const pair<double, Candidate>& b) {
                    return a.first < b.first;
                });
                sort(notSuitableCandidates.begin(), notSuitableCandidates.end(), [](const pair<double, Candidate>& a, const pair<double, Candidate>& b) {
                    return a.first < b.first;
                });
                
                cout << "\nCandidates suitable for job \"" << targetJob.title << "\" (from most to least suitable):" << endl;
                if (suitableCandidates.empty()) {
                    cout << "No suitable candidate found!" << endl;
                } else {
                    for (const auto &entry : suitableCandidates) {
                        cout << "Candidate: " << entry.second.name 
                             << ", ID: " << entry.second.id 
                             << ", Hybrid Cost: " << entry.first << endl;
                    }
                }
                
                cout << "\nCandidates not suitable for job \"" << targetJob.title << "\" (Hybrid Cost above threshold):" << endl;
                if (notSuitableCandidates.empty()) {
                    cout << "All candidates are suitable for this job." << endl;
                } else {
                    for (const auto &entry : notSuitableCandidates) {
                        cout << "Candidate: " << entry.second.name 
                             << ", ID: " << entry.second.id 
                             << ", Hybrid Cost: " << entry.first << endl;
                    }
                }
            }
        }else if (choice == 6) {
            string candName, jobTitle;
            int candId, jobId;
            
            // Get candidate details.
            cout << "Enter Candidate name: ";
            cin.ignore();
            getline(cin, candName);
            cout << "Enter Candidate ID: ";
            cin >> candId;
            
            // Get job details.
            cin.ignore();
            cout << "Enter Job title: ";
            getline(cin, jobTitle);
            cout << "Enter Job ID: ";
            cin >> jobId;
            
            // Find candidate.
            bool candidateFound = false;
            Candidate targetCandidate;
            for (const auto &cand : candidates) {
                if (cand.name == candName && cand.id == candId) {
                    targetCandidate = cand;
                    candidateFound = true;
                    break;
                }
            }
            
            // Find job.
            bool jobFound = false;
            Job targetJob;
            for (const auto &job : jobs) {
                if (toLower(job.title) == toLower(jobTitle) && job.id == jobId) {
                    targetJob = job;
                    jobFound = true;
                    break;
                }
            }
            
            if (!candidateFound) {
                cout << "Candidate not found!" << endl;
            } else if (!jobFound) {
                cout << "Job not found!" << endl;
            } else {
                // Generate match.
                MatchResult result = generateMatch(targetCandidate, targetJob);
                
                // Display match details and suitability analysis.
                displayMatchAnalysis(targetCandidate, targetJob, result.cost, result.suitable);
            }
        }else if (choice == 7) {
            // Option 7: Search for candidates by multiple skills.
            searchCandidatesByMultipleSkills(candidates);
        }else if(choice==8){
            cout << "Exiting the system..." << endl;
            break;
        } else {
            cout << "Invalid choice. Please enter a valid option." << endl;
        }
    }

    return 0;
}