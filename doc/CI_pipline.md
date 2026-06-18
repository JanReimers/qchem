### 1. GitHub Actions Self-Hosted Runner (Recommended if code is on GitHub)
If your code is already on GitHub, this is the easiest integration. You keep the GitHub UI but run the jobs on your hardware.

*   **Setup Steps**:
    1.  Go to your GitHub Repo **Settings** > **Actions** > **Runners** > **New self-hosted runner**.
    2.  Download the runner package on your Linux box and extract it.
    3.  Configure it with the provided token:
        ```bash
        ./config.sh --url https://github.com/USER/REPO --token YOUR_TOKEN
        ```
    4.  Install it as a service: `sudo ./svc.sh install` and `sudo ./svc.sh start`.
*   **Workflow**: Create a `.github/workflows/ci.yml` file in your repo:
    ```yaml
    jobs:
      build:
        runs-on: self-hosted  # Targets your Linux box
        steps:
          - uses: actions/checkout@v4
          - run: mkdir build && cd build && cmake .. && make && ctest
    ```
*   **Benefit**: No new UI to learn; your existing GitHub workflow triggers the build on your machine automatically.

### 2. Jenkins (Standalone Automation Server)
Jenkins is a dedicated Java-based server ideal if you want a standalone dashboard independent of GitHub/GitLab.

*   **Installation**:
    ```bash
    # Install Java
    sudo apt install openjdk-17-jre -y
    # Add Jenkins repo and install
    sudo wget -O /usr/share/keyrings/jenkins-keyring.asc https://pkg.jenkins.io/debian-stable/jenkins.io-2023.key
    echo "deb [signed-by=/usr/share/keyrings/jenkins-keyring.asc] https://pkg.jenkins.io/debian-stable binary/" | sudo tee /etc/apt/sources.list.d/jenkins.list
    sudo apt update && sudo apt install jenkins -y
    sudo systemctl start jenkins
    ```
*   **Configuration**:
    1.  Access `http://localhost:8080` and complete the setup wizard.
    2.  Install the **"CMake"** and **"JUnit"** (for gtest XML output) plugins.
    3.  Create a **Pipeline** job that points to your repo.
    4.  Define your pipeline in a `Jenkinsfile`:
        ```groovy
        pipeline {
            agent any
            stages {
                stage('Build & Test') {
                    steps {
                        sh 'mkdir build && cd build && cmake .. && make'
                        sh 'cd build && ctest --output-on-failure'
                    }
                }
            }
        }
        ```

### 3. GitLab CI (Self-Hosted)
If you host your code on a self-hosted GitLab instance, you use **GitLab Runners**.

*   **Installation**:
    ```bash
    sudo apt install gitlab-runner
    sudo gitlab-runner register # Follow prompts to link to your project
    sudo gitlab-runner start
    ```
*   **Configuration**: Create a `.gitlab-ci.yml` in your repo root:
    ```yaml
    stages:
      - test
    cpp_tests:
      stage: test
      script:
        - mkdir build && cd build
        - cmake ..
        - make
        - ctest --output-on-failure
    ```

### Critical Step: CTest Integration
Regardless of the tool, ensure your CMake configuration outputs results in a format the CI can parse (e.g., JUnit XML for Jenkins/GitHub):

```cmake
# In CMakeLists.txt
enable_testing()
# Add this to generate XML reports for CI dashboards
set(CTEST_OUTPUT_ON_FAILURE TRUE)
```
Then run tests in your pipeline with:
```bash
ctest --output-on-failure -T Test --no-compress-output
```
This generates XML reports in the `Testing` directory, which Jenkins or GitHub Actions can archive to show visual test trends.


