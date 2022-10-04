
$(function () {
    $('#toggle > button').click(function () {
        var ix = $(this).index();

        $('#left').toggle(ix === 0);
        $('#right').toggle(ix === 1);
    });

});

var desktopBtn = $("#desktop");
var mobileBtn = $("#mobile");
var body = $('body');

desktopBtn.on('click', function () {
    body.addClass('large-screen');
    togglePrimaryButtonStyle($(this));
})

mobileBtn.on('click', function () {
    body.removeClass('large-screen');
    togglePrimaryButtonStyle($(this));
})

function togglePrimaryButtonStyle(el) {
    var sibling = el.parent('.btn-group').siblings('.btn-group').find('.btn');
    el.addClass('btn-primary');
    sibling.removeClass('btn-primary').addClass('btn-default');
}