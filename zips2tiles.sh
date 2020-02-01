#!/usr/bin/env bash

ls | grep UTM33 | grep 133 | parallel ./norgezip2tiles.sh ./t2.py {}
