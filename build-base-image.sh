#!/bin/bash
# Script to build and push the base image with RevBayes compiled with bundled Boost
# This only needs to be run when RevBayes needs to be updated (rare)
# NOTE: Requires lib/revbayes/ submodule to be initialized and updated first

set -e

# Cleanup on failure
trap 'echo "Build failed! Cleaning up..."; docker rmi "${IMAGE_NAME}" 2>/dev/null || true' ERR

# Pre-flight checks
echo "Performing pre-flight checks..."

# Check if Docker is running
if ! docker info >/dev/null 2>&1; then
    echo "❌ Error: Docker is not running or not installed"
    echo "Please start Docker Desktop or install Docker first"
    exit 1
fi

# Check available disk space (need at least 10GB)
AVAILABLE_SPACE=$(df -k . | tail -1 | awk '{print $4}')
REQUIRED_SPACE=$((10 * 1024 * 1024))  # 10GB in KB
if [ "$AVAILABLE_SPACE" -lt "$REQUIRED_SPACE" ]; then
    echo "⚠️  Warning: Low disk space detected"
    echo "Available: $((AVAILABLE_SPACE / 1024 / 1024))GB, Recommended: 10GB+"
    echo "Continue anyway? (y/N)"
    read -r response
    if [[ ! "$response" =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo "✅ Pre-flight checks passed"
echo ""

# Generate date-based tag
DATE_TAG=$(date +%Y-%m-%d)
IMAGE_NAME="quay.io/matsengrp/linearham:${DATE_TAG}-base-image"

echo "=========================================="
echo "Building base image: ${IMAGE_NAME}"
echo "=========================================="
echo ""
echo "This will take 20-40 minutes as it compiles RevBayes from source with bundled Boost."
echo "Resource limits: 8GB memory"
echo "Press Ctrl+C to cancel, or wait 5 seconds to continue..."
sleep 5

# Build the base image with resource limits
echo ""
echo "Building base image..."
docker build --platform linux/amd64 \
  --memory=8g \
  --progress=plain \
  -f Dockerfile.base \
  -t "${IMAGE_NAME}" .

echo ""
echo "=========================================="
echo "Build complete!"
echo "=========================================="
echo ""
echo "Base image built: ${IMAGE_NAME}"
echo ""
echo "Next steps:"
echo "1. Test the image:"
echo "   docker run --rm ${IMAGE_NAME} rb --version"
echo ""
echo "2. Push to quay.io:"
echo "   docker push ${IMAGE_NAME}"
echo ""
echo "3. Update the main Dockerfile FROM line to:"
echo "   FROM ${IMAGE_NAME}"
echo ""
echo "4. Commit and push the Dockerfile change"
echo ""
