import pytest
from pathlib import Path
import subprocess

# Define the root of the repository and the directory containing Python scripts
repo_root = Path(__file__).resolve().parent.parent
scripts_directory = repo_root / "saltOptimization"

@pytest.fixture
def python_scripts():
    # Using .rglob to recursively find all .py files
    return list(scripts_directory.rglob('*.py'))

def test_python_scripts(python_scripts):
    for script_path in python_scripts:
        # Get script path relative to repository root
        relative_script_path = script_path.relative_to(repo_root)
        print(f"Testing script: {relative_script_path}")

        with subprocess.Popen(['python', str(script_path)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
            try:
                # Read the output line by line in real-time
                for line in proc.stdout:
                    print(line, end='')  # end='' to avoid double newlines

                # Wait for the subprocess to finish
                proc.wait()

                # Check if the process exited with a non-zero status code
                if proc.returncode != 0:
                    pytest.fail(f"Script at {relative_script_path} failed with return code: {proc.returncode}")

            except Exception as e:
                pytest.fail(f"An error occurred while running script {relative_script_path}: {e}")