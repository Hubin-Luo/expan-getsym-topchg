#!/bin/bash

#get the DMI energy per site
qstate=$1
engfile=$2
natom=20
lino=`awk '$1=="iq=" && $2=='$qstate'{print NR}' $engfile`
awk 'BEGIN{n=0};
    {
     if(NR=='$lino'+2+3*n && n<'$natom')
       {
	edm=0
	for(i=1;i<=12;i++)
           {
	    edm+=$i
           }
	printf "%4d %12.4f\n", n+1, edm
        n+=1
       }
    }' $engfile > edm-$qstate.dat
