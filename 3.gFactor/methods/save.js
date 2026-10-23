// Shared by every task. On the Sussex server data goes to that task's PHP script as before. Anywhere else
// (e.g. GitHub Pages) it goes to DataPipe, which stores each file in the study's Zenodo draft record;
// ../import_datapipe/import_datapipe.py later runs the same PHP checks on the downloaded files.
const DATAPIPE_EXPERIMENT_ID = 'kF37MjUSvYs9' //from the DataPipe experiment dashboard (Zenodo)
const ON_SUSSEX = window.location.hostname === 'users.sussex.ac.uk'
const SAVE_TIMEOUT = 15000 //ms per attempt - DataPipe replies quickly, queueing the upload itself if storage is slow
const SAVE_RETRY_DELAYS = [2000, 5000] //ms before each automatic retry of a failed attempt
const RESEARCHER_EMAIL = 'm.lovell@sussex.ac.uk'

// Only moves on (next_task) once the data is confirmed saved. Connection/server problems are retried, then the
// participant gets a "Try again" button and a download of their responses to email, so nothing is lost silently.
function postData(php_url, file_name, exp_data, next_task){
    let url, body
    if(ON_SUSSEX){
        url = php_url
        body = JSON.stringify({file_name: file_name, exp_data: exp_data})
    } else {
        url = 'https://pipe.jspsych.org/api/data/' //https://pipe.jspsych.org/docs/api
        body = JSON.stringify({experimentID: DATAPIPE_EXPERIMENT_ID, filename: file_name + '.json', data: JSON.stringify(exp_data)})
    }
    let attempt = 0
    function send(){
        attempt++
        const xhr = new XMLHttpRequest()
        xhr.onloadend = function() { //fires after success, error, abort or timeout
            const outcome = saveOutcome(xhr)
            if(outcome === 'saved'){
                hideSaveProblem()
                next_task()
            } else if(outcome === 'retry' && attempt <= SAVE_RETRY_DELAYS.length){
                showSaveProblem('trying', file_name, exp_data, send)
                setTimeout(send, SAVE_RETRY_DELAYS[attempt-1])
            } else {
                showSaveProblem(outcome, file_name, exp_data, send)
            }
        }
        xhr.open('POST', url, true)
        xhr.timeout = SAVE_TIMEOUT
        xhr.setRequestHeader('Content-Type', 'application/json')
        xhr.send(body)
    }
    send()
}

// 'saved', 'retry' (connection or server problem - worth trying again) or 'stuck' (trying again won't help)
function saveOutcome(xhr){
    if(xhr.status >= 200 && xhr.status < 300){ return 'saved' } //PHP replies 200; DataPipe 201 stored or 202 queued
    if(xhr.status === 400){
        let error = ''
        try { error = JSON.parse(xhr.responseText).error } catch(e) {}
        if(error === 'FILE_EXISTS'){ return 'saved' } //an earlier attempt already got through
        return 'stuck' //e.g. INVALID_DATA or SESSION_LIMIT_REACHED
    }
    return 'retry' //status 0 (no connection or timed out) or a 5xx server error
}

function showSaveProblem(state, file_name, exp_data, retry){
    let box = document.getElementById('save_problem')
    if(!box){
        box = document.createElement('div')
        box.id = 'save_problem'
        box.style.cssText = 'position:fixed;top:0;left:0;right:0;bottom:0;z-index:9999;display:flex;align-items:center;' +
            'justify-content:center;padding:16px;background:rgba(255,255,255,0.97);color:#000;font:18px/1.5 sans-serif;text-align:center'
        document.body.appendChild(box)
    }
    if(state === 'trying'){
        box.innerHTML = '<p>Saving your responses is taking longer than usual - trying again...</p>'
        return
    }
    const why = state === 'stuck' ? 'Your responses could not be saved.'
                                  : 'Your responses could not be saved. This is usually a problem with the internet connection.'
    box.innerHTML = '<div style="max-width:34em"><p>' + why + ' Please keep this page open.</p>' +
        '<p><button type="button" id="save_retry" style="font-size:18px;padding:8px 20px">Try again</button></p>' +
        '<p>If it still doesn\'t work, <a id="save_download" href="#">download your responses</a> and email the file to ' +
        RESEARCHER_EMAIL + ' so your participation still counts.</p></div>'
    document.getElementById('save_retry').onclick = function() {
        box.innerHTML = '<p>Saving your responses...</p>'
        retry()
    }
    const link = document.getElementById('save_download')
    link.href = URL.createObjectURL(new Blob([JSON.stringify(exp_data)], {type: 'application/json'}))
    link.download = file_name + '.json' //same name and contents DataPipe would store, so import_datapipe.py can check it
}

function hideSaveProblem(){
    const box = document.getElementById('save_problem')
    if(box){ box.remove() }
}
