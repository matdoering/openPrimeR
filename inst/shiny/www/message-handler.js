Shiny.addCustomMessageHandler("testmessage",
    function(message) {
        alert(JSON.stringify(message));
    }
); 
Shiny.addCustomMessageHandler("jsCode",
    function(message) {
        eval(message.value);
    }
);
