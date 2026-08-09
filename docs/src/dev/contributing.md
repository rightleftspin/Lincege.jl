## Contributing

Additions to Lincege are welcome in the following forms:

1. **Issues** – Report bugs or suggest new features.
2. **Pull Requests** – Submit code changes or documentation updates.

Please ensure all your edits pass the local tests located in `test/runtests.jl` before submitting.

### Pull Request Process

We use the feature branch workflow for contributions. Please follow these steps:

1. **Fork the Repository**
   Fork the Lincege repository to your own GitHub account.

2. **Clone Locally**
   Clone your fork to your computer and navigate into the directory:
   ```bash
   git clone https://github.com
   cd Lincege.jl
   ```

3. **Create a Feature Branch**
   Create a new branch for your specific changes. Use a descriptive name:
   ```bash
   git checkout -b feature/your-feature-name
   ```

4. **Make Changes and Test**
   Implement your changes. Ensure everything works by running the tests:
   ```julia
   using Pkg; Pkg.test("Lincege")
   ```
   *(Alternatively, run `julia test/runtests.jl` from your terminal)*

5. **Commit and Push**
   Commit your changes with clear messages and push the branch to your fork:
   ```bash
   git add .
   git commit -m "Add a brief description of your changes"
   git push origin feature/your-feature-name
   ```

6. **Open a Pull Request**
   Go to the original Lincege repository on GitHub. Click **New Pull Request** and select your feature branch to submit.
