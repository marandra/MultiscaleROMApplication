#!/bin/bash

for SCRIPT in tmp_*.bash
do
  pueue add -- bash $SCRIPT
  sleep 0.05
done
