Qualtrics.SurveyEngine.addOnload(function () {
    //https://api.qualtrics.com/api-reference/ZG9jOjg0MDczOA-api-reference
    //document.addEventListener('touchstart',function(e){console.log(e)})
    //document.addEventListener('touchmove',function(e){console.log(e)})
    
    const qthis = this;
    qthis.hideNextButton();
    const resize = document.getElementById("resize")
    document.body.insertBefore(resize,document.body.firstChild)

    const resizer = document.getElementById("resizer")
    const box = document.getElementById("box")
    const resize_button = document.getElementById("resize_button")

    box.style.width = box.clientWidth
    box.style.height = box.clientHeight
    resizer.addEventListener('mousedown', registerResizer, false);
    resizer.addEventListener('touchstart', handleResize, false);
    resize_button.addEventListener('click', resizeButton, false);

    function handleResize(e){
        e.preventDefault();
        if(e.type==='touchstart'){ 
            e = e.changedTouches[0] 
            resizer.addEventListener('touchmove', handleResize, false);
        } else if(e.type==='touchmove'){
            e = e.changedTouches[0]
        }
        const box_location = box.getBoundingClientRect();
        box.style.width = (e.pageX - box_location.x) + 'px';
        box.style.height = (e.pageY - box_location.y) + 'px';
        if(e.pageX >= window.innerWidth){ box.style.width = (window.innerWidth-1)+'px'}
        if(e.pageY >= window.innerHeight){ box.style.height = (window.innerHeight-1)+'px' }
    }

    function registerResizer(e){ //mouse
        window.addEventListener('mousemove', handleResize, false);
        window.addEventListener('mouseup', unregisterResizer, false);
    }

    function unregisterResizer(e){ //remove mouse
        window.removeEventListener('mousemove', handleResize, false);
        window.removeEventListener('mouseup', unregisterResizer, false);
    }

    function resizeButton(e){ //finish resizer and forward
        const credit_card = [8.56,5.398] //https://myaurochs.com/blogs/news/whats-a-credit-cards-size
        const box = document.getElementById("box")
        const px_cm_w = box.clientWidth/credit_card[0]
        const px_cm_h = box.clientHeight/credit_card[1]
        const px_cm = (px_cm_w+px_cm_h)/2
        Qualtrics.SurveyEngine.setEmbeddedData('px_cm', px_cm) //must be set in survey flow https://stackoverflow.com/questions/35732280/how-to-set-embedded-data-in-qualtrics-with-javascript
        resize.remove();
        qthis.clickNextButton(); //a few backups to push to next qualtrics screen just incase.
        qthis.showNextButton();
        jQuery("#NextButton").click();
        return
    }

})

Qualtrics.SurveyEngine.addOnReady(function(){});
Qualtrics.SurveyEngine.addOnUnload(function(){});