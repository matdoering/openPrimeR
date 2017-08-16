/*
 $(function () {
   $(document).on('shown.bs.tooltip', function (e) {
      setTimeout(function () {
        $(e.target).tooltip('destroy');
      }, 10);
   });
});
*/

shinyBS.addTooltip = function(id, type, opts) {
  var $id = shinyBS.getTooltipTarget(id);
  /* increase the overall tooltip delay */
  var dopts = {html: true, delay: { "show": 1000, "hide": 200 }};
  opts = $.extend(opts, dopts);
    //alert(opts);
  if(type == "tooltip") {
    $id.tooltip("destroy");
    $id.tooltip(opts);
  } else if(type == "popover") {
    $id.popover("destroy");
    $id.popover(opts);
  }
  
}
