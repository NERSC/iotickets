<?php

//$user = "rharred_ro";
//$password = "Eyef409c2b0";
//$host = "webdbr";
//$database = "jobs";

$user = "rharred";
$password = "journey";
$host = "c3po";
$database = "queuestats";

mysql_connect($host,$user,$password);
@mysql_select_db($database) or die( "Unable to select database");
//$query="SELECT AVG(start - queued) AS MYVAR FROM summary where (hostname = 'edison' AND completion > '1404082489' )";
//$query="SELECT * FROM jobQueue where (NodeSize = $myNodes AND Walltime = '$myWallTime')";
$query="SELECT * FROM jobQueue";
$result=mysql_query($query);
$num=mysql_numrows($result);
mysql_close();

//echo $num;

echo "[";

$imfirst = TRUE;

$i=0;
while ($i < $num) {

  $QPosition=mysql_result($result,$i,"QPosition");
  $StartTime=mysql_result($result,$i,"StartTime");
  if ($QPosition != 0 && !empty($StartTime)) {
    if ($imfirst) {
      $imfirst = FALSE;
      echo "{";
    } else {
      echo ",{";
    }
    echo "\"QPosition\" : $QPosition,";
    echo "\"StartTime\" : \"$StartTime\",";
    $TimeCalled=mysql_result($result,$i,"TimeCalled");
    echo "\"TimeCalled\" : \"$TimeCalled\",";
    $Backlog=mysql_result($result,$i,"Backlog");
    echo "\"Backlog\" : $Backlog,";
    $Rsv=mysql_result($result,$i,"Rsv");
    echo "\"Rsv\" : \"$Rsv\",";
    $Priority=mysql_result($result,$i,"Priority");
    echo "\"Priority\" : \"$Priority\",";
    $JobID=mysql_result($result,$i,"JobID");
    echo "\"JobID\" : $JobID,";
    $Historical=mysql_result($result,$i,"Historical");
    echo "\"Historical\" : \"$Historical\"";
    echo "}";
  }
  $i++;
}

echo "]";

?> 
