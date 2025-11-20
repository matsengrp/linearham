#!/bin/bash
# Quick local test script for Dockerfile changes
# This mimics what CI does without running the full workflow

set -e

echo "Building Docker image locally..."
docker build --platform linux/amd64 -t linearham-test:local .

echo ""
echo "✅ Docker build succeeded!"
echo ""
echo "To run tests locally:"
echo "  docker run linearham-test:local sh -c '/linearham/_build/test/test'"
echo "  docker run linearham-test:local sh -c '/linearham/test.sh'"
echo ""
echo "To clean up:"
echo "  docker rmi linearham-test:local"
