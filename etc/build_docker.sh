#! /bin/bash

set -e

IMAGE_NAME=pcoc
DOCKERFILE_DIR=.
REPO=carinerey/$IMAGE_NAME
docker build -t $REPO $DOCKERFILE_DIR

if [[ $1 == "push_yes" ]]
then 
    docker push $REPO
fi
