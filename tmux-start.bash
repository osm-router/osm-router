#!/bin/sh
cd ./src/
SESSION="osm-router"

tmux -2 new-session -d -s $SESSION

tmux new-window -t $SESSION:1 -k -n main
tmux send-keys -t $SESSION:1 'vim Utils.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe ../dists.cfg' C-m
tmux send-keys -t $SESSION:1 '1gt' C-m
tmux split-window -h
tmux send-keys -t $SESSION:1 'vim Router.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe Router.c++' C-m
tmux select-pane -t 0

cd ..
tmux new-window -t $SESSION:2 -n makefile
tmux select-window -t $SESSION:2
tmux send-keys -t $SESSION:2 'vim .travis.yml' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe tmux-start.bash' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe CMakeLists.txt' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe README.md' C-m
tmux split-window -h
tmux select-pane -t 0

cd ./extract_profiles/
tmux new-window -t $SESSION:3 -n extract_profile

tmux select-window -t $SESSION:3
tmux send-keys -t $SESSION:3 'vim extract.c++' C-m
tmux split-window -h
tmux send-keys -t $SESSION:3 'vim Makefile' C-m
tmux split-window -v
tmux select-pane -t 0

tmux attach -t $SESSION
