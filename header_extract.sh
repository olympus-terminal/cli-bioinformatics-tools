#!/bin/bash

fgrep \> >> "$LINE".heads

sed -i 's/>//g' "$LINE".heads
