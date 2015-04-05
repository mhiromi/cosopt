#!/bin/bash
for host in `awk '{print $1;}' ~/mpd.hosts`
do
	rsh $host killall python2.3
	rsh $host killall mpi_cosopt
	rsh $host rm /tmp/mpd2.console_`whoami`
done
