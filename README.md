
# Candidate-Job-Matching

## Overview
This project implements an automated job matching system designed as a case study for Design and Analysis of Algorithms. The system efficiently pairs job candidates with suitable job positions by integrating advanced algorithmic techniques such as the Hungarian Algorithm for optimal assignment and string matching methods (Jaccard Similarity with tokenization and Smith-Waterman) for evaluating textual compatibility.

## Problem Definition
The objective of the project is to develop a system that can automatically match candidates to job openings based on multiple criteria. Instead of relying solely on keyword searches, this system evaluates candidate profiles and job descriptions using:
- **Structured criteria** (skills, certifications, experience, salary expectations, location, soft skills)
- **Unstructured text similarity** (resume and job description matching)

This approach ensures a globally optimal assignment that minimizes the overall mismatch cost.

## Methodology

### 1. Network Flow Matching using the Hungarian Algorithm
- **Assignment Problem Formulation:**  
  A cost matrix is constructed where each cell represents the mismatch cost between a candidate and a job, computed by a hybrid cost function.
- **Hungarian Algorithm:**  
  The algorithm iteratively reduces the matrix (row and column reductions, covering zeros, and adjustments) to find the optimal one-to-one assignment. It guarantees a globally optimal solution in O(n³) time.
- **Handling Unequal Sets:**  
  Dummy rows/columns (with a very high cost) are added to balance the matrix when the number of candidates and jobs differ.

### 2. String Matching Algorithms
- **Jaccard Similarity with Tokenization:**  
  - **Purpose:** Compares candidate skills, certifications, and soft skills with job requirements by converting text to lowercase tokens and filtering out common stop words.

- **Smith-Waterman Algorithm:**  
  - **Purpose:** Computes local sequence alignment between a candidate’s resume and the job description, yielding a normalized similarity score that is then translated into a cost.

### 3. Hybrid Cost Function
- **Combination:**  
  The overall cost for a candidate-job pair is computed as a weighted sum of the structured cost (skills, experience, salary, certifications, soft skills) and the text-based cost (resume vs. job description similarity).
- **Minimization:**  
  Lower cost values indicate better matches. The Hungarian Algorithm minimizes the total cost to achieve an optimal assignment.
  
This combination ensures a balance between accuracy and efficiency for typical candidate-job matching problems.

## Data Format
The system reads candidate and job information from a JSON file. Here is an example structure:

### Candidate Structure
```json
{
  "id": 1,
  "name": "Alice",
  "resume": "Experienced software engineer with expertise in C++ and Java.",
  "skills": ["C++", "Java", "SQL"],
  "cgpa": 3.8,
  "experience": 5,
  "certificates": [
    {
      "title": "Certified C++ Developer",
      "issuer": "Microsoft"
    },
    {
      "title": "Java Professional Certification",
      "issuer": "Oracle"
    }
  ],
  "expectedSalary": 80000,
  "preferredLocation": "NYC",
  "softSkills": ["communication", "teamwork", "problem-solving"]
}
```

### Job Structure
```json
{
  "id": 101,
  "title": "Software Engineer",
  "description": "Looking for an experienced engineer skilled in C++ and SQL.",
  "requiredSkills": ["C++", "SQL"],
  "min_cgpa": 3.5,
  "min_experience": 4,
  "requiredCertifications": ["Certified C++ Developer"],
  "offeredSalary": 80000,
  "location": "NYC",
  "requiredSoftSkills": ["communication", "teamwork"]
}
```

## Installation & Usage

### Prerequisites
- **C++ Compiler** (g++ recommended)
- **JSON for Modern C++ Library** (nlohmann/json)

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/akshayks13/Candidate-Job-Matching.git
   ```
   
2. Ensure the `json.hpp` file from the [nlohmann/json](https://github.com/nlohmann/json) library is in the `include/` directory or accessible in your include path.

### Compilation
```bash
g++ -o app app.cpp -std=c++11
```

### Running the Application
```bash
./app
```
The program will display a menu to list candidates, jobs, generate match reports, search by skills, and more.

## Project Workflow
1. **Data Parsing:** Reads and parses candidate and job data from JSON.
2. **Cost Matrix Construction:** Computes a hybrid cost for every candidate-job pair.
3. **Matrix Padding:** Balances the cost matrix for unequal numbers of candidates and jobs.
4. **Assignment:** Uses the Hungarian Algorithm to find the optimal assignment.
5. **Match Analysis:** Provides detailed analysis on why a candidate is suitable (or not) for a job.
 
