#!/bin/bash
# create_symlink.sh

# Replace 'some_directory' with the actual directory you want to link to
THERMOCHIMICA="$HOME/thermochimica"
THERMOCHIMICA_LINK="thermochimica"
MSTDB_TP="$HOME/mstdb-tp"
MSTDB_TP_LINK="mstdb-tp"

# Create the symlink
ln -sfn "$THERMOCHIMICA" "$THERMOCHIMICA_LINK"
ln -sfn "$MSTDB_TP" "$MSTDB_TP_LINK"