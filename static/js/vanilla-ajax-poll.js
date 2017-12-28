
// The polling function
function poll(fn, timeout, interval) {
    var endTime = Number(new Date()) + (timeout || 2000);
    interval = interval || 100;

    var checkCondition = function(resolve, reject) { 
        var ajax = fn();
        // dive into the ajax promise
        ajax.then( function(response){
            // If the condition is met, we're done!
            if(response.data.var == true) {
                resolve(response.data.var);
            }
            // If the condition isn't met but the timeout hasn't elapsed, go again
            else if (Number(new Date()) < endTime) {
                setTimeout(checkCondition, interval, resolve, reject);
            }
            // Didn't match and too much time, reject!
            else {
                reject(new Error('timed out for ' + fn + ': ' + arguments));
            }
        });
    };

    return new Promise(checkCondition);
}

// Usage: get something via ajax
//poll(function() {
//	return axios.get('something.json');
//}, 2000, 150).then(function() {
//    // Polling done, now do something else!
//}).catch(function() {
//    // Polling timed out, handle the error!
//});

//The polling function
function poll(fn, timeout, interval) {
    var endTime = Number(new Date()) + (timeout || 2000);
    interval = interval || 100;

    var checkCondition = function(resolve, reject) {
        // If the condition is met, we're done! 
        var result = fn();
        if(result) {
            resolve(result);
        }
        // If the condition isn't met but the timeout hasn't elapsed, go again
        else if (Number(new Date()) < endTime) {
            setTimeout(checkCondition, interval, resolve, reject);
        }
        // Didn't match and too much time, reject!
        else {
            reject(new Error('timed out for ' + fn + ': ' + arguments));
        }
    };

    return new Promise(checkCondition);
}

//Usage:  ensure element is visible
poll(function() {
	return document.getElementById('lightbox').offsetWidth > 0;
}, 2000, 150).then(function() {
    // Polling done, now do something else!
}).catch(function() {
    // Polling timed out, handle the error!
});


