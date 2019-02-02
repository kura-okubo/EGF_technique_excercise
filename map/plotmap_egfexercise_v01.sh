#!/usr/local/bin/zsh
psbasemap -JM12c -R-130/-115/30/50 -B5NseW -X5 -Y4 -K -P> map.eps

pscoast -J -R -W -K -O -Gtan -S -I1 >> map.eps
pscoast -J -R -W -K -O -Bg5a10f5   -Lf-120/32.5/32.5/400 >> map.eps
psxy station.dat -J -R -Si0.5 -G0/205/0 -W -K -O >> map.eps
psmeca focaldata.dat -Jm -R -Sa1 -G255/0/0 -L -P -O >>map.eps


#psbasemap -JM12c -R-135/-120/30/41  -B5NseW -X5 -Y4 -K -P> inset.eps
#pscoast -J -R -W -K -O -Bg1a2f1 >> inset.eps
#psmeca focaldata.dat -Jm -R -Sa1 -G255/0/0 -O >>inset.eps
