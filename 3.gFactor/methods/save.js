// Shared by every task. On the Sussex server data goes to that task's PHP script as before. Anywhere else
// (e.g. GitHub Pages) it goes to DataPipe, which stores each file in Google Drive; ../import_datapipe/import_datapipe.py
// later runs the same PHP checks on those files and uploads them to the Sussex server.
const DATAPIPE_EXPERIMENT_ID = 'REPLACE_WITH_DATAPIPE_EXPERIMENT_ID' //from the DataPipe experiment dashboard
const ON_SUSSEX = window.location.hostname === 'users.sussex.ac.uk'

function postData(php_url, file_name, exp_data, next_task){
    const xhr = new XMLHttpRequest()
    xhr.onloadend = function() { next_task() } //move on to next task - onloadend also fires on error/timeout so participants can't get stuck
    let body
    if(ON_SUSSEX){
        xhr.open('POST', php_url, true)
        body = JSON.stringify({file_name: file_name, exp_data: exp_data})
    } else {
        xhr.open('POST', 'https://pipe.jspsych.org/api/data/', true) //https://pipe.jspsych.org/docs/api
        body = JSON.stringify({experimentID: DATAPIPE_EXPERIMENT_ID, filename: file_name + '.json', data: JSON.stringify(exp_data)})
    }
    xhr.timeout = 30000 //ms - stops a hung request blocking onloadend; DataPipe allows up to 60 s while it writes to Drive
    xhr.setRequestHeader('Content-Type', 'application/json')
    xhr.send(body)
}
