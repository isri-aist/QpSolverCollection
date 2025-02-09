#!/bin/sh
set -e

apt-get update
apt-mark hold ros-*
apt-get upgrade -y

basic_dep="git \
           curl \
           nano \
           vim \
           python3-pip \
           python3-colcon-common-extensions \
           wget \
           ninja-build \
           x11-apps"

apt-get install -y $basic_dep
apt-get update
apt-get upgrade -y

mkdir /root/catkin_ws
cd /root/catkin_ws

# Install eigen-qld
# curl -1sLf \
#   'https://dl.cloudsmith.io/public/mc-rtc/stable/setup.deb.sh' \
#   | sudo -E bash

# apt install libeigen-qld-dev python3-eigen-qld python3-eigen-qld -y

# Clean up
apt-get autoremove -y
apt-get clean -y