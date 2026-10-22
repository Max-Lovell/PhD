<?php
// Runs one of the study's own save scripts (e.g. methods/breath/breath.php) from the command line, exactly as the
// Sussex server would for a browser POST, except the request body comes from stdin and the file is written to
// <out_dir> instead of ../../gfactor/tasks. Called by import_datapipe.py:
//   php run_php.php <task.php> <out_dir> < body.json
// The task script exits silently on anything it rejects, so "no file written" means "rejected".
error_reporting(E_ALL & ~E_DEPRECATED & ~E_WARNING); // FILTER_SANITIZE_STRING is deprecated in PHP 8.1+; missing fields warn before the script exits

[, $__php_file, $__out_dir] = $argv;
$__src = file_get_contents($__php_file);
$__replacements = [ // each must appear exactly once, so an edited PHP file fails loudly instead of writing somewhere unexpected
    "file_get_contents('php://input')" => 'stream_get_contents(STDIN)',
    '"../../gfactor/tasks/$file_name.json"' => '$__out_dir."/$file_name.json"',
];
foreach ($__replacements as $__from => $__to) {
    if (substr_count($__src, $__from) !== 1) {
        fwrite(STDERR, "Expected exactly one $__from in $__php_file\n");
        exit(2);
    }
    $__src = str_replace($__from, $__to, $__src);
}
// span.php has a live debug dump of $GLOBALS; on the command line that includes every environment variable
$__src = str_replace('file_put_contents("globals.log", print_r($GLOBALS,true));', '', $__src);

$_SERVER['REQUEST_METHOD'] = 'POST';
$_SERVER['HTTP_ORIGIN'] = 'https://users.sussex.ac.uk';
$_SERVER['HTTP_USER_AGENT'] = getenv('GFACTOR_USER_AGENT') ?: ''; // survey.php stores this as user_agent
eval('?>' . $__src);
