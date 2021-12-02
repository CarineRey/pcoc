#! /bin/bash

set -e

IMAGE_NAME=pcoc
TAG=v1.1.0
DOCKERFILE_DIR=.
REPO=carinerey/$IMAGE_NAME
docker build -t $REPO:$TAG -f etc/Dockerfile $DOCKERFILE_DIR
docker build -t $REPO -f etc/Dockerfile $DOCKERFILE_DIR
if [[ $1 == "push_yes" ]]
then 
    docker push $REPO:$TAG
    docker push $REPO
fi
