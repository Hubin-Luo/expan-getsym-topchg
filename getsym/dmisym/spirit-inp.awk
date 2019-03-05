#!/usr/bin/awk -f

#convert the output to spirit input

BEGIN{
      ox=0.0; oy=0.0; oz=0.0
      v1=0.0; v2=0.0; v3=0.0
      p1=0.0; p2=0.0; p3=0.0
     }
     {
      if(NF==5){
        ox=$3; oy=$4; oz=$5
      }
      else{
        v1=$3; v2=$4; v3=$5
        p1=ox+v1; p2=oy+v2; p3=oz+v3
        if(p1<0){p1=int(p1)-1}
	if(p1>=0){p1=int(p1)}
        if(p2<0){p2=int(p2)-1}
        if(p2>=0){p2=int(p2)}
        if(p3<0){p3=int(p3)-1}
        if(p3>=0){p3=int(p3)}
	printf "%4d %6d %6d %4d %4d %16.9f %16.9f %16.9f %16.9f\n", \
        $1-1, $2-1, p1, p2, p3, $6, $7, $8, $9*13.6057
      }
     }
