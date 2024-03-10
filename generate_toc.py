import os

# -------------------------------------------------------------------------------
# This is a script for automatically updating the table of contents of the main
# README file in the LEUplus directory, please don't update this unless changing
# the TOC is your intention.
# -------------------------------------------------------------------------------

def generate_toc(directory, level=0):
    toc = ""
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        if os.path.isdir(item_path):
            readme_path = os.path.join(item_path, "README.md")
            if os.path.isfile(readme_path):
                indent = "  " * level
                toc += f"{indent}- [{item}]({os.path.relpath(readme_path, '.')})\n"
                toc += generate_toc(item_path, level + 1)
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

# Specify the path to your main README file
readme_path = "README.md"

# Generate the Table of Contents
toc = generate_toc(".")

# Update the Table of Contents in the README file if there are changes
update_toc_in_readme(readme_path, toc)