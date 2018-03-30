#!/bin/bash

# Add local user
# Either use the LOCAL_USER_ID if passed in at runtime or
# fallback

USER_ID=${LOCAL_USER_ID:-9001}
Xvfb :1 -screen 0 1024x768x16 &
export DISPLAY=:1

if [ -n "$CWD" ]
then
W_DIR=$CWD
elif [ -n "$SHARED_DIR" ]
then
W_DIR=$SHARED_DIR
fi

W_DIR=${W_DIR:-/data}

echo "Working directory : $W_DIR"
cd $W_DIR

if [ -n "$LOCAL_USER_ID" ] 
then
USER_ID=${LOCAL_USER_ID:-9001}
echo "Starting with UID : $USER_ID"
useradd --shell /bin/bash -u $USER_ID -o -c "" -g sudo -m user
export HOME=/home/user
echo "user ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
exec /usr/sbin/gosu user "$@"
else
echo "Starting with UID : root"
exec "$@"
fi
