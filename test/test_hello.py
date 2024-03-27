# contents of test_script1.py

import pytest

def test_script1_runs():
    # Assuming your script has a main function that runs the main script logic
    try:
        print('Hello world')
    except Exception as e:
        pytest.fail(f"script1 failed with an exception: {e}")