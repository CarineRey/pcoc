#! /bin/bash

set -e

IMAGE_NAME=pcoc
#TAG=master
DOCKERFILE_DIR=.
REPO=carinerey/$IMAGE_NAME:$TAG
docker build -t $REPO -f etc/Dockerfile $DOCKERFILE_DIR

if [[ $1 == "push_yes" ]]
then 
    docker push $REPO
fi
