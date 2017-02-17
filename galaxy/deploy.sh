#!/bin/bash
# this script uploads galaxy tools in the current directory to the galaxy toolsheds
set -e

CURRENT_DIR=$(dirname $0)

echo "Installing planemo ... again"
pip install planemo

echo "Deploying to Testtoolshed ..."
planemo shed_update -r --force_repository_creation -t testtoolshed --shed_key_from_env TTS_KEY "$CURRENT_DIR"

echo "Deploying to Toolshed ..."
planemo shed_update -r --force_repository_creation -t toolshed --shed_key_from_env TS_KEY $CURRENT_DIR

echo "Successfully deployed to toolsheds"
exit 0
