node=0
offset=0
filename=ranks
dims=$1

echo "# Manual mapping of MPI ranks to CPU cores ($dims-D)" > $filename
echo "#" >> $filename

if [ $dims -eq 2 ]
then
echo "# Pure MPI: 2D level 4 SGCT with 2^13x2^13 grid points with two layers, -p 64, -q 32" >> $filename
echo "# proposed process layout with grid_0 = 16x4, grid_1 = 8x8, grid_2 = 8x8, grid_3 = 4x16," >> $filename
echo "# grid_4 = 8x4, grid_5 = 8x4, grid_6 = 4x8, full = 16x16" >> $filename
elif [ $dims -eq 3 ]
then
echo "# Pure MPI: 3D level 4 SGCT with 2^8x2^8x2^8 grid points with three layers, -p 64, -q 32, -r 16" >> $filename
echo "# proposed process layout with grid_0 = 16x2x2, grid_1 = 8x4x2, grid_2 = 4x8x2, grid_3 = 2x16x2," >> $filename
echo "# grid_4 = 8x2x4, grid_5 = 4x4x4, grid_6 = 2x8x4, grid_7 = 4x2x8, grid_8 = 2x4x8, grid_9 = 2x2x16," >> $filename
echo "# grid_10 = 8x2x2, grid_11 = 4x4x2, grid_12 = 2x8x2, grid_13 = 4x2x4, grid_14 = 2x4x4," >> $filename 
echo "# grid_15 = 2x2x8, grid_16 = 4x2x2, grid_17 = 2x4x2, grid_18 = 2x2x4, full = 8x8x8" >> $filename
fi

echo "# ======================================================================================" >> $filename
echo "#" >> $filename

while read line;do
        ##################################################################
        ## 2D                                                           ##
        ##################################################################
        if [ $dims -eq 2 ]
        then
	   if [ $node -eq 0 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=0
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 4 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 5 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 6 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 7 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 20 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 21 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 22 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 23 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 1 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=8
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 4 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 5 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 6 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 7 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 20 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 21 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 22 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 23 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 2 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=32
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 4 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 5 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 6 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 7 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 20 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 21 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 22 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 23 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 3 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=40
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 4 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 5 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 6 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 7 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 20 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 21 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 22 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 23 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 4 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=64
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 5 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=68
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 6 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=96
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 7 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=100
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 8 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=128
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 9 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=132
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 10 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=160
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 11 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=164
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 12 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=192
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 21 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 29 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 13 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=194
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 21 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 29 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 14 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=224
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 21 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 29 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 15 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=226
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 21 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 29 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 16 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=256
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 17 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=260
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 18 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=288
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 19 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=292
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 20 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=320
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 21 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 29 + $offset`=$nodename slot=1:7" >> $filename
	   elif [ $node -eq 21 ]
	   then
	   nodename=$(echo "$line" | awk '{print $1;}')
	   offset=322
	   echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
	   echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
	   echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
	   echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
	   echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
	   echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
	   echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
	   echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
	   echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
	   echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
	   echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
	   echo "rank `expr 21 + $offset`=$nodename slot=1:3" >> $filename
	   echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
	   echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
	   echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
	   echo "rank `expr 29 + $offset`=$nodename slot=1:7" >> $filename
	   fi
        ##################################################################
        ## 3D                                                           ##
        ##################################################################
        elif [ $dims -eq 3 ]
        then
           if [ $node -eq 0 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=0
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 7 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 35 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 37 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 38 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 39 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 1 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=8
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 7 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 35 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 37 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 38 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 39 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 2 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=16
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 7 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 35 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 37 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 38 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 39 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 3 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=24
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 7 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 35 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 37 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 38 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 39 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 4 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=64
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 35 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 41 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 42 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 43 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 5 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=68
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 35 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 41 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 42 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 43 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 6 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=80
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 35 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 41 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 42 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 43 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 7 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=84
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 35 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 41 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 42 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 43 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 8 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=128
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 37 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 41 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 44 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 45 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 9 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=130
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 37 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 41 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 44 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 45 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 10 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=144
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 37 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 41 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 44 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 45 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 11 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=146
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 37 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 41 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 44 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 45 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 12 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=192
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 38 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 42 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 44 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 46 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 13 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=193
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 38 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 42 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 44 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 46 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 14 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=208
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 38 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 42 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 44 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 46 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 15 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=209
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 38 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 42 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 44 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 46 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 16 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=256
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 19 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 35 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 49 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 50 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 51 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 17 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=260
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 19 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 35 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 49 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 50 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 51 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 18 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=264
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 19 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 35 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 49 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 50 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 51 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 19 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=268
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 19 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 35 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 49 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 50 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 51 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 20 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=320
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 21 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 37 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 49 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 52 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 53 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 21 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=322
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 21 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 37 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 49 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 52 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 53 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 22 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=328
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 21 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 37 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 49 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 52 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 53 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 23 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=330
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 21 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 37 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 49 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 52 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 53 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 24 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=384
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 22 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 38 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 50 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 52 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 54 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 25 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=385
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 22 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 38 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 50 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 52 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 54 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 26 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=392
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 22 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 38 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 50 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 52 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 54 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 27 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=393
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 22 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 38 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 50 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 52 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 54 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 28 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=448
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 25 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 41 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 49 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 56 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 57 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 29 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=450
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 25 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 41 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 49 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 56 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 57 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 30 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=452
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 25 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 41 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 49 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 56 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 57 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 31 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=454
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 25 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 33 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 41 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 49 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 56 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 57 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 32 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=512
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 26 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 42 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 50 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 56 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 58 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 33 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=513
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 26 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 42 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 50 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 56 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 58 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 34 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=516
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 26 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 42 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 50 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 56 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 58 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 35 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=517
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 26 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 34 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 42 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 50 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 56 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 58 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 36 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=576
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 44 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 52 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 56 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 60 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 37 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=577
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 44 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 52 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 56 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 60 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 38 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=578
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 44 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 52 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 56 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 60 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 39 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=579
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 32 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 36 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 40 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 44 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 48 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 52 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 56 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 60 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 40 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=640
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 41 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=644
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 11 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 19 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 26 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 27 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 42 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=672
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 21 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 29 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 43 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=674
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 21 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 29 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 44 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=704
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 22 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 26 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 30 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 45 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=705
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 22 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 26 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 30 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 46 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=736
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 21 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 29 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 47 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=738
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 13 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 17 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 21 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 25 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 29 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 48 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=768
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 22 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 26 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 30 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 49 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=769
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 22 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 26 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 30 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 50 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=800
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 22 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 26 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 30 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 51 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=801
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 16 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 18 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 20 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 22 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 24 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 26 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 28 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 30 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 52 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=832
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 7 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 11 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 13 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 15 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 53 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=848
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 7 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 11 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 13 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 15 + $offset`=$nodename slot=1:7" >> $filename
           elif [ $node -eq 54 ]
           then
           nodename=$(echo "$line" | awk '{print $1;}')
           offset=864
           echo "rank `expr 0 + $offset`=$nodename slot=0:0" >> $filename
           echo "rank `expr 1 + $offset`=$nodename slot=0:1" >> $filename
           echo "rank `expr 2 + $offset`=$nodename slot=0:2" >> $filename
           echo "rank `expr 3 + $offset`=$nodename slot=0:3" >> $filename
           echo "rank `expr 4 + $offset`=$nodename slot=0:4" >> $filename
           echo "rank `expr 5 + $offset`=$nodename slot=0:5" >> $filename
           echo "rank `expr 6 + $offset`=$nodename slot=0:6" >> $filename
           echo "rank `expr 7 + $offset`=$nodename slot=0:7" >> $filename
           echo "rank `expr 8 + $offset`=$nodename slot=1:0" >> $filename
           echo "rank `expr 9 + $offset`=$nodename slot=1:1" >> $filename
           echo "rank `expr 10 + $offset`=$nodename slot=1:2" >> $filename
           echo "rank `expr 11 + $offset`=$nodename slot=1:3" >> $filename
           echo "rank `expr 12 + $offset`=$nodename slot=1:4" >> $filename
           echo "rank `expr 13 + $offset`=$nodename slot=1:5" >> $filename
           echo "rank `expr 14 + $offset`=$nodename slot=1:6" >> $filename
           echo "rank `expr 15 + $offset`=$nodename slot=1:7" >> $filename
           fi
        fi
        ((node++))
done < hostfile
# Format of hostfile is
#
# hostname0 slots=16
# hostname1 slots=16

chmod 755 $filename

