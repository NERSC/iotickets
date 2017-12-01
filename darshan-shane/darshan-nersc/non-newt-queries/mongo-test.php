<?php require 'mongodb-creds.php';

echo json_encode(array('connected' => $mongo->connected));

?>
