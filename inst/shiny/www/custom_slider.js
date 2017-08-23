$(document).ready(function() {
    /**
    Custom slider labels
    */
    // Convert the region position (e.g. from -60 to 40) to leader/exon position
    function positionToLabel(pos) {
        label = "";
        if (pos < 0) {
            label = pos; 
        } else if (pos == 0) {
            label = "target region start";
        } else {
            label = "+" + pos;
        }
        return label;
    }
    function percentageFormat(ratio) {
        ratio = ratio * 100;
        ratio = ratio + "%"
        return ratio;
    }
    var slider_customize_allowed_fw = $("#individual_allowed_regions_fw").ionRangeSlider({
        prettify: positionToLabel,
        force_edges: true,
        grid: true,
        grid_num: 4,
    });
    var slider_customize_allowed_rev = $("#individual_allowed_regions_rev").ionRangeSlider({
        prettify: positionToLabel,
        force_edges: true,
        grid: true,
        grid_num: 4
    });
    var slider_gc_ratio_allowed = $("#allowed_gc_ratio").ionRangeSlider({
        prettify: percentageFormat,
        force_edges: true,
        grid: true,
        grid_num: 4
    });
    var slider_gc_ratio_limit = $("#limit_allowed_gc_ratio").ionRangeSlider({
        prettify: percentageFormat,
        force_edges: true,
        grid: true,
        grid_num: 4
    });
    var slider_primer_efficiency_allowed = $("#allowed_primer_efficiency").ionRangeSlider({
        prettify: percentageFormat,
        force_edges: true,
        grid: true,
        grid_num: 4
    });
    var slider_primer_specificity_allowed = $("#allowed_primer_specificity").ionRangeSlider({
        prettify: percentageFormat,
        force_edges: true,
        grid: true,
        grid_num: 4
    });
    var slider_primer_specificity_limit = $("#limit_allowed_primer_specificity").ionRangeSlider({
        prettify: percentageFormat,
        force_edges: true,
        grid: true,
        grid_num: 4
    });
    var slider_allowed_other_binding_ratio = $("#allowed_other_binding_ratio").ionRangeSlider({
        prettify: percentageFormat,
        force_edges: true,
        grid: true,
        grid_num: 4
    });

    var slider_required_conservation = $("#required_conservation").ionRangeSlider({
        prettify: percentageFormat,
        force_edges: true,
        grid: true,
        grid_num: 4
    });

    var slider_required_opti_cvg = $("#required_opti_cvg").ionRangeSlider({
        prettify: percentageFormat,
        force_edges: true,
        grid: true,
        grid_num: 4
    });
    var slider_allowed_model_FPR = $("#allowed_coverage_model").ionRangeSlider({
        prettify: percentageFormat,
        force_edges: true,
        grid: true,
        grid_num: 4
    });

}) 


