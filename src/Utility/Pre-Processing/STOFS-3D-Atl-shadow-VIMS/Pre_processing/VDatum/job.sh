#/bin/bash
java -jar vdatum.jar ihorz:NAD83_2011 ivert:navd88:m:height ohorz:igs14:geo:deg overt:xgeoid20b:m:height -file:txt:space,1,2,3,skip0:hgrid_stofs3d_ches_del.txt:result region:5 &
java -jar vdatum.jar ihorz:NAD83_2011 ivert:navd88:m:height ohorz:igs14:geo:deg overt:xgeoid20b:m:height -file:txt:space,1,2,3,skip0:hgrid_stofs3d_inland_1.txt:result region:4 &
java -jar vdatum.jar ihorz:NAD83_2011 ivert:navd88:m:height ohorz:igs14:geo:deg overt:xgeoid20b:m:height -file:txt:space,1,2,3,skip0:hgrid_stofs3d_inland_2.txt:result region:4 &
java -jar vdatum.jar ihorz:NAD83_2011 ivert:navd88:m:height ohorz:igs14:geo:deg overt:xgeoid20b:m:height -file:txt:space,1,2,3,skip0:hgrid_stofs3d_inland_3.txt:result region:4 &
java -jar vdatum.jar ihorz:NAD83_2011 ivert:navd88:m:height ohorz:igs14:geo:deg overt:xgeoid20b:m:height -file:txt:space,1,2,3,skip0:hgrid_stofs3d_inland_4.txt:result region:4 &
