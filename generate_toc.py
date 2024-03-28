import os, sys
import re

# Get the script name from sys.argv
script_name = sys.argv[0]

# Get the absolute path of the script
script_dir = Path(os.path.abspath(script_name)).resolve().parent

# -------------------------------------------------------------------------------
# This is a script for automatically updating the table of contents of the main
# README file in this repository, please don't update this unless changing
# the TOC is your intention.
# -------------------------------------------------------------------------------

def load_gitignore_patterns(gitignore_path=str(script_dir / '.gitignore')):
    patterns = []
    with open(gitignore_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith('#'):
                # Convert the .gitignore pattern to a regular expression
                pattern = re.escape(line).replace(r'\*', '.*').replace(r'\?', '.')
                patterns.append(re.compile(pattern))
    return patterns

def is_ignored(path, ignore_patterns):
    for pattern in ignore_patterns:
        if pattern.fullmatch(path):
            return True
    return False

def generate_toc(directory, ignore_patterns, level=0):
    toc = ""
    for item in sorted(os.listdir(directory)):
        item_path = os.path.join(directory, item)
        if os.path.isdir(item_path) and not is_ignored(item, ignore_patterns):
            readme_path = os.path.join(item_path, "README.md")
            if os.path.isfile(readme_path):
                indent = "  " * level
                toc += f"{indent}- [{item}]({os.path.relpath(readme_path, '.')})\n"
                toc += generate_toc(item_path, ignore_patterns, level + 1)
    return toc

def update_toc_in_readme(readme_path, new_toc):
    with open(readme_path, "r") as file:
        content = file.read()

    toc_start = content.find("## Table of Contents")
    if toc_start != -1:
        toc_end = content.find("##", toc_start + 1)
        if toc_end == -1:
            toc_end = len(content)
        existing_toc = content[toc_start:toc_end].strip()
        if existing_toc == f"## Table of Contents\n\n{new_toc}":
            print("Table of Contents is up to date. No changes made.")
            return
        updated_content = content[:toc_start] + "## Table of Contents\n\n" + new_toc.strip() + "\n\n" + content[toc_end:]
        with open(readme_path, "w") as file:
            file.write(updated_content)
        print("Table of Contents updated successfully.")
    else:
        print("Error: Table of Contents section not found in the README file.")

# Load .gitignore patterns
gitignore_patterns = load_gitignore_patterns()

# Specify the path to your main README file
readme_path = "README.md"

# Generate the Table of Contents
toc = generate_toc(".", gitignore_patterns)

# Update the Table of Contents in the README file if there are changes
update_toc_in_readme(readme_path, toc)