<?php include "data.php";

//$user = $_GET['username'];

if( isset($_GET['oldest']) && isset($_GET['newest']) ) {

	$cursor = queryList($_GET['username'], $_GET['oldest'], $_GET['newest']);

} else if ( isset($_GET['jobid'])) {

	$cursor = queryListForJobID($_GET['username'], $_GET['jobid']);

} else {

	$cursor = queryList($_GET['username']);

}

$sorted = array();

foreach($cursor as $document) {

	$job = array(

		'_id' => $document['_id'],
		'jobid' => $document['jobid'],
		'end_time' => $document['end_time'],
		'host' => $document['host'],
		'exe' => $document['exe']

	);

	$sorted[] = $job;

}

$length = sizeof($sorted);

for($i = 0; $i < $length; $i++) {

	for($j = $i; $j < $length; $j++) {

		if( $sorted[$j]['end_time'] > $sorted[$i]['end_time'] ) {

			$temp = $sorted[$i];
			$sorted[$i] = $sorted[$j];
			$sorted[$j] = $temp;

		}

	}

}

echo json_encode($sorted);

// courtesy of Shane Pearlman

?>


