Usage of DOUG's web service wrapper:

Serverside:
1) Have Java installed
2) Install Tomcat
3) enable WebDAV in Tomcat
4) edit Tomcat user file, so that Authentication is necessary.
   Give authenticated users right to read and write in WebDAV directory.
5) download and deploy Axis (not Axis2) webapp in Tomcat
6) edit constants in ee.ut.math.doug.DougService as necessary
6) compile serverside code of DOUG web service wrapper and put it into 
   $CATALINA_HOME/webapps/axis/WEB-INF/classes.
7) use deploy.wsdd on axis webapp to register web serive.   
8) compile DOUG for server platform and put it into WebDAV directory.

Clientside:
1) put data files into WebDAV directory.
2) compile client side code of DOUG web service wrapper.
3) run 'java DougWSClient [-h hostname | -p port]'
