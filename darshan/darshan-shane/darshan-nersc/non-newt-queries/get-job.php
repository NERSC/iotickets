<?php include "data.php";

$job = $_GET['_id'];

$cursor = query($job);

echo json_encode($cursor->getNext());

?>


