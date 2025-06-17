#!/bin/bash

# clone repositories
git clone https://github.com/andrewlhicks/mymesh.git
git clone https://github.com/andrewlhicks/sympyplus.git
git clone https://github.com/andrewlhicks/q-tensor-3d.git

# install packages
pip install -e mymesh
pip install -e sympyplus
pip install -e q-tensor-3d

# get current directory
CURRENT_DIR=$(pwd)

NEW_PATH=$CURRENT_DIR/q-tensor-3d/scripts

if ! echo "$PATH" | grep -q ":$NEW_PATH:"; then
    # add scripts to PATH for current shell session
    PATH="$PATH:$NEW_PATH"
    echo "Added $NEW_PATH to PATH"

    # add scripts to PATH in "default" shell
    if [[ $SHELL = /bin/bash ]]; then
        CONFIG_FILE=$HOME/.bashrc
    elif [[ $SHELL = /bin/zsh ]]; then
        CONFIG_FILE=$HOME/.zshrc
    else
        echo "Cannot export PATH because bash or zsh was not detected. Please add $NEW_PATH to your PATH manually."
    fi

    echo "Updating $CONFIG_FILE to include $NEW_PATH in PATH"
    echo "export PATH=\"\$PATH:$NEW_PATH\"" >> $CONFIG_FILE
else
    echo "$NEW_PATH is already in PATH"
fi

# create petsc configuration file
echo "-options_left 0" > $HOME/.petscrc