#!/bin/sh
if [ $# -lt 1 ]; then
    option=1
else
    option=$1
fi

cd ./src/
SESSION="osm-router"

tmux -2 new-session -d -s $SESSION

tmux new-window -t $SESSION:1 -k -n main
if [ $option -eq 1 ]
then
    tmux send-keys -t $SESSION:1 'vim Utils.h' C-m
    tmux send-keys -t $SESSION:1 ':' 'tabe ../profile.cfg' C-m
    tmux send-keys -t $SESSION:1 '1gt' C-m
else
    tmux send-keys -t $SESSION:1 'vim xml-parser.h' C-m
    tmux send-keys -t $SESSION:1 ':' 'tabe xml-parser.c++' C-m
    tmux send-keys -t $SESSION:1 ':' 'tabe Graph.h' C-m
    tmux send-keys -t $SESSION:1 ':' 'tabe Graph.c++' C-m
fi
tmux split-window -h
if [ $option -eq 1 ]
then
    tmux send-keys -t $SESSION:1 'vim Router.h' C-m
    tmux send-keys -t $SESSION:1 ':' 'tabe Router.c++' C-m
    tmux select-pane -t 0
else
    tmux send-keys -t $SESSION:1 'vim Router-test.h' C-m
    tmux send-keys -t $SESSION:1 ':' 'tabe Router-test.c++' C-m
fi

cd ..
tmux new-window -t $SESSION:2 -n makefile
tmux select-window -t $SESSION:2
tmux send-keys -t $SESSION:2 'vim .travis.yml' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe tmux-start.bash' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe CMakeLists.txt' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe README.md' C-m
tmux split-window -h
if [ $option -ne 1 ]
then
    tmux send-keys -t $SESSION:2 'git st' C-m
    tmux send-keys -t $SESSION:2 'cd build/' C-m
fi
tmux select-pane -t 0

if [ $option -eq 1 ]
then
    cd ./extract_profiles/
    tmux new-window -t $SESSION:3 -n extract_profile

    tmux select-window -t $SESSION:3
    tmux send-keys -t $SESSION:3 'vim extract.c++' C-m
    tmux split-window -h
    tmux send-keys -t $SESSION:3 'vim Makefile' C-m
    tmux split-window -v
    tmux select-pane -t 0
fi

tmux select-window -t $SESSION:1

tmux attach -t $SESSION
