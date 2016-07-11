{\rtf1\ansi\ansicpg1252\cocoartf1343\cocoasubrtf160
{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural

\f0\fs36 \cf0 README
\fs24 \
\
Run:\
++++\
\
python ValueIteration.py\
\
This will complete both value and policy iteration, output the run times, and then ask you to query a state (you can select from policy or value iteration). A state is defined by the number of open rooms, and by the total number of people at the start of the hour (=people in the waiting room + arrival of new people).\
\
For help run:\
++++++++++\
\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural
\cf0 python ValueIteration.py -h\
\
This will describe all the different parameters.\
\
Example\
+++++++\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural

\f1\fs22 \cf0 \CocoaLigature0 python ValueIteration.py -N 10 -M 5 --E1=10 --E2=2 --E3=80 -W 5 -P 300 --dist=1 -l 3 -h}