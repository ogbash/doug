-- MySQL dump 10.9
--
-- Host: localhost    Database: dougtests
-- ------------------------------------------------------
-- Server version	4.1.13

--
-- Current Database: `dougtests`
--

CREATE DATABASE /*!32312 IF NOT EXISTS*/ `dougtests` /*!40100 DEFAULT CHARACTER SET latin1 */;

USE `dougtests`;

--
-- Table structure for table `testresultfiles`
--

DROP TABLE IF EXISTS `testresultfiles`;
CREATE TABLE `testresultfiles` (
  `ID` int(11) NOT NULL auto_increment,
  `testresult_ID` int(11) default NULL,
  `name` varchar(200) default NULL,
  `content` mediumblob,
  PRIMARY KEY  (`ID`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--
-- Table structure for table `testresults`
--

DROP TABLE IF EXISTS `testresults`;
CREATE TABLE `testresults` (
  `ID` int(11) NOT NULL auto_increment,
  `testrun_ID` int(11) default NULL,
  `starttime` datetime default NULL,
  `endtime` datetime default NULL,
  `method` int(11) default NULL,
  `solver` int(11) default NULL,
  `status` int(11) default NULL,
  `nproc` int(11) default NULL,
  `errortext` text,
  `name` varchar(60) default NULL,
  inputtype int,
  levels int,
  executable varchar(16),
  PRIMARY KEY  (`ID`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--
-- Table structure for table `testruns`
--

DROP TABLE IF EXISTS `testruns`;
CREATE TABLE `testruns` (
  `ID` int(11) NOT NULL auto_increment,
  `servername` varchar(64) default NULL,
  `starttime` datetime default NULL,
  `endtime` datetime default NULL,
  svnrevision int,
  fcompiler varchar(32),
  mpi varchar(32),
  errortext text,
  PRIMARY KEY  (`ID`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

