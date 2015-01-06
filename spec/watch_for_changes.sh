#!/bin/bash
DIRECTORY_TO_OBSERVE=`pwd`
echo 'watching' $DIRECTORY_TO_OBSERVE
function block_for_change {
  inotifywait -r \
    -e modify,move,create,delete \
    $DIRECTORY_TO_OBSERVE
}

echo 'running' $1
$1

while block_for_change; do
  echo 'running' $1
  $1
done
