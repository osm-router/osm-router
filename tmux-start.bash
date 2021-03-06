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
    tmux send-keys -t $SESSION:1 'vim Router.h' C-m
    tmux send-keys -t $SESSION:1 '1gt' C-m
else
    tmux send-keys -t $SESSION:1 'vim xml-parser.h' C-m
    tmux send-keys -t $SESSION:1 ':' 'tabe xml-parser.cpp' C-m
    tmux send-keys -t $SESSION:1 ':' 'tabe Graph.h' C-m
    tmux send-keys -t $SESSION:1 ':' 'tabe Graph.cpp' C-m
fi
tmux split-window -h
if [ $option -eq 1 ]
then
    tmux send-keys -t $SESSION:1 'vim Router.cpp' C-m
    tmux select-pane -t 0
else
    tmux send-keys -t $SESSION:1 'vim stochastic-sp.h' C-m
    tmux send-keys -t $SESSION:1 ':' 'tabe stochastic-sp.cpp' C-m
    tmux send-keys -t $SESSION:1 ':' 'tabe stochastic-test.h' C-m
    tmux send-keys -t $SESSION:1 ':' 'tabe stochastic-test.cpp' C-m
    tmux send-keys -t $SESSION:1 '1gt' C-m
fi

cd ..
tmux new-window -t $SESSION:2 -n makefile
tmux select-window -t $SESSION:2
tmux send-keys -t $SESSION:2 'vim .travis.yml' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe tmux-start.bash' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe CMakeLists.txt' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe README.md' C-m
tmux split-window -h
if [ $option -eq 1 ]
then
    tmux select-pane -t 0
else
    tmux send-keys -t $SESSION:2 'git st' C-m
    tmux send-keys -t $SESSION:2 'cd build/' C-m
fi

tmux select-window -t $SESSION:1

tmux attach -t $SESSION
