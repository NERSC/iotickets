<?php include "data.php";

$user = $_GET['username'];

if( isset($_GET['oldest']) && isset($_GET['newest']) ) {

  $cursor = queryAggregate($user, $_GET['oldest'], $_GET['newest']);

} else if ( isset($_GET['jobid'])) {

  $cursor = queryAggregateForJob($user, $_GET['jobid']);

} else {

  $cursor = queryAggregate($user);

}

$numOfDocs = $cursor->count(true);

$aggData = array(

  "_id" => 'n/a',
  "start_time" => 'n/a',
  "end_time" => 'n/a',
  "tot_time" => 0,
  "io_time" => 0,
  "nprocs" => 'n/a',
  "host" => 'n/a',
  "exe" => 'n/a',
  "uid" => 'n/a',
  "total_CP_BYTES_READ" => 0,
  "total_CP_BYTES_WRITTEN" => 0,
  "total_CP_SIZE_READ_0_100" => 0,
  "total_CP_SIZE_READ_100_1K" => 0,
  "total_CP_SIZE_READ_1K_10K" => 0,
  "total_CP_SIZE_READ_10K_100K" => 0,
  "total_CP_SIZE_READ_100K_1M" => 0,
  "total_CP_SIZE_READ_1M_4M" => 0,
  "total_CP_SIZE_READ_4M_10M" => 0,
  "total_CP_SIZE_READ_10M_100M" => 0,
  "total_CP_SIZE_READ_100M_1G" => 0,
  "total_CP_SIZE_READ_1G_PLUS" => 0,
  "total_CP_SIZE_WRITE_0_100" => 0,
  "total_CP_SIZE_WRITE_100_1K" => 0,
  "total_CP_SIZE_WRITE_1K_10K" => 0,
  "total_CP_SIZE_WRITE_10K_100K" => 0,
  "total_CP_SIZE_WRITE_100K_1M" => 0,
  "total_CP_SIZE_WRITE_1M_4M" => 0,
  "total_CP_SIZE_WRITE_4M_10M" => 0,
  "total_CP_SIZE_WRITE_10M_100M" => 0,
  "total_CP_SIZE_WRITE_100M_1G" => 0,
  "total_CP_SIZE_WRITE_1G_PLUS" => 0,
  "total_CP_F_POSIX_READ_TIME" => 0,
  "total_CP_F_POSIX_WRITE_TIME" => 0,
  "total_CP_F_POSIX_META_TIME" => 0,
  "agg_perf_by_cumul" => 0,
  "agg_perf_by_open" => 0,
  "agg_perf_by_open_lastio" => 0,
  "agg_perf_by_slowest" => 0,
  "total_bytes" => 0


);

$keys = array(
  "total_CP_BYTES_READ",
  "total_CP_BYTES_WRITTEN",
  "total_CP_SIZE_READ_0_100",
  "total_CP_SIZE_READ_100_1K",
  "total_CP_SIZE_READ_1K_10K",
  "total_CP_SIZE_READ_10K_100K",
  "total_CP_SIZE_READ_100K_1M",
  "total_CP_SIZE_READ_1M_4M",
  "total_CP_SIZE_READ_4M_10M",
  "total_CP_SIZE_READ_10M_100M",
  "total_CP_SIZE_READ_100M_1G",
  "total_CP_SIZE_READ_1G_PLUS",
  "total_CP_SIZE_WRITE_0_100",
  "total_CP_SIZE_WRITE_100_1K",
  "total_CP_SIZE_WRITE_1K_10K",
  "total_CP_SIZE_WRITE_10K_100K",
  "total_CP_SIZE_WRITE_100K_1M",
  "total_CP_SIZE_WRITE_1M_4M",
  "total_CP_SIZE_WRITE_4M_10M",
  "total_CP_SIZE_WRITE_10M_100M",
  "total_CP_SIZE_WRITE_100M_1G",
  "total_CP_SIZE_WRITE_1G_PLUS",
  "total_CP_F_POSIX_READ_TIME",
  "total_CP_F_POSIX_WRITE_TIME",
  "total_CP_F_POSIX_META_TIME",
  "agg_perf_by_cumul",
  "agg_perf_by_open",
  "agg_perf_by_open_lastio",
  "agg_perf_by_slowest",
  "total_bytes"
);

foreach($cursor as $document) {

	foreach($keys as $key) {

		if(is_numeric($document[$key]))
			$aggData[$key] += $document[$key];

	}
        $aggData['tot_time'] += $document['end_time'] - $document['start_time'];
        $aggData['io_time'] += $document['total_bytes']/1024/1024/$document['agg_perf_by_slowest'];

}

foreach($keys as $key) {

	$aggData[$key] = floor($aggData[$key]/$numOfDocs*1000)/1000;

}

echo json_encode($aggData);

?>
