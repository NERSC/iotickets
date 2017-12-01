<?php require 'mongodb-creds.php';

$collection = $mongo->darshan->slog;

$fields = array(
	"jobid" => true,
	"start_time" => true,
	"end_time" => true,
	"nprocs" => true,
	"host"  => true,
	"exe"  => true,
	"uid"  => true,
	"total_CP_BYTES_READ" => true,
	"total_CP_BYTES_WRITTEN" => true,
	"total_CP_SIZE_READ_0_100" => true,
	"total_CP_SIZE_READ_100_1K" => true,
	"total_CP_SIZE_READ_1K_10K" => true,
	"total_CP_SIZE_READ_10K_100K" => true,
	"total_CP_SIZE_READ_100K_1M" => true,
	"total_CP_SIZE_READ_1M_4M" => true,
	"total_CP_SIZE_READ_4M_10M" => true,
	"total_CP_SIZE_READ_10M_100M" => true,
	"total_CP_SIZE_READ_100M_1G" => true,
	"total_CP_SIZE_READ_1G_PLUS" => true,
	"total_CP_SIZE_WRITE_0_100" => true,
	"total_CP_SIZE_WRITE_100_1K" => true,
	"total_CP_SIZE_WRITE_1K_10K" => true,
	"total_CP_SIZE_WRITE_10K_100K" => true,
	"total_CP_SIZE_WRITE_100K_1M" => true,
	"total_CP_SIZE_WRITE_1M_4M" => true,
	"total_CP_SIZE_WRITE_4M_10M" => true,
	"total_CP_SIZE_WRITE_10M_100M" => true,
	"total_CP_SIZE_WRITE_100M_1G" => true,
	"total_CP_SIZE_WRITE_1G_PLUS" => true,
	"total_CP_F_POSIX_READ_TIME" => true,
	"total_CP_F_POSIX_WRITE_TIME" => true,
	"total_CP_F_POSIX_META_TIME" => true,
	"agg_perf_by_cumul"  => true,
	"agg_perf_by_open"  => true,
	"agg_perf_by_open_lastio"  => true,
	"agg_perf_by_slowest"  => true,
	"total_bytes" => true
);

function queryList($user, $oldest, $newest) {

	global $collection;

	$regexObj = new MongoRegex("/^$user\_/");

	$query = array();

	$query = array_merge($query, array('_id' => $regexObj));

	if( isset($oldest) && isset($newest) ) {

		$query = array_merge( $query, array(

			'$and' => array(

				array(

					'end_time' => array('$gte' => (int) $oldest)

				),

				array( 

					'end_time' => array('$lte' => (int) $newest)

				)

			)

		));

	}

	$fields = array(

		'jobid' => true,
		'end_time' => true,
		'host' => true,
		'exe' => true

	);

	$cursor = $collection->find($query,$fields);

	return $cursor;
	
}

function queryListForJobID($user, $jobid) {

	global $collection;

	$regexObj = new MongoRegex("/^$user\_/");

	$query = array();

	$query = array_merge($query, array('_id' => $regexObj));

	$query = array_merge( $query, array( 'jobid' => (int) $jobid ));

	$fields = array(

		'jobid' => true,
		'end_time' => true,
		'host' => true,
		'exe' => true

	);

	$cursor = $collection->find($query,$fields);

	return $cursor;
	
}

function query($job) {

	global $collection;

	$query = array(

		'_id' => $job

	);

	global $fields;

	$cursor = $collection->find($query,$fields);

	return $cursor;

}

function queryTimeRange($user) {

	global $collection;

	$regexObj = new MongoRegex("/^$user\_/");

	$query = array();

	$query = array_merge($query, array('_id' => $regexObj));
        $query = array_merge($query, array('end_time' => array('$exists' => true)));

	$fields = array(
		'_id' => false,
		'end_time' => true
	);

	$oldest = $collection->find($query,$fields)->sort(array(
			'end_time' => 1
		))->limit(1);

	$oldest = iterator_to_array($oldest);

	if(sizeof($oldest) > 0) {
		$oldest = $oldest[0]['end_time'];
	} else {
		$oldest = null;
	}

	$newest = $collection->find($query,$fields)->sort(array(
			'end_time' => -1
		))->limit(1);

	$newest = iterator_to_array($newest);

	if(sizeof($newest) > 0) {
		$newest = $newest[0]['end_time'];
	} else {
		$newest = null;
	}

	$range = array(
		'oldest' => $oldest,
		'newest' => $newest
	);

	return $range;

}

function queryAggregate($user, $oldest, $newest) {

	global $collection;

	$regexObj = new MongoRegex("/^$user\_/");

	$query = array();

	$query = array_merge($query, array('_id' => $regexObj));

	if( isset($oldest) && isset($newest) ) {

		$query = array_merge( $query, array(

			'$and' => array(

				array(

					'end_time' => array('$gte' => (int) $oldest)

				),

				array( 

					'end_time' => array('$lte' => (int) $newest)

				)

			)

		));

	}

	global $fields;

	$cursor = $collection->find($query,$fields);

	return $cursor;

}

function queryAggregateForJob($user, $jobid) {

	global $collection;

	$regexObj = new MongoRegex("/^$user\_/");

	$query = array();

	$query = array_merge($query, array('_id' => $regexObj));

	$query = array_merge($query, array( 'jobid' => (int) $jobid ));

	global $fields;

	$cursor = $collection->find($query,$fields);

	return $cursor;

}

?>
