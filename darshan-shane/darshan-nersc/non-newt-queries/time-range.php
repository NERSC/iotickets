<?php include "data.php";

$user = $_GET['username'];

$range = queryTimeRange($user);

echo json_encode($range);

?>
