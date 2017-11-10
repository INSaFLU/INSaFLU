;(function ( $, window, document, undefined ) {

	/*
		DjangoModal HandlerResponse
	*/
	var DjangoModalHandlerResponse = function(){}

	DjangoModalHandlerResponse.prototype.handle = function(response){
		var res = null;
		if(response.type === 'normal'){
			res = this.handleNormalResponse(response);
		} else {
			res = this.handleRedirectResponse(response);
		}
		return res
	}

	DjangoModalHandlerResponse.prototype.handleNormalResponse = function(response){
		return response.content;
	}

	DjangoModalHandlerResponse.prototype.handleRedirectResponse = function(response){
		window.location = response.redirect_to;
	}

	DjangoModalHandlerResponse.prototype.handleErrorResponse = function(jqXHR, textStatus, errorThrown){

	}


	/*
		DjangoModal AjaxForm
	
	*/

	var DjangoModalAjaxForm = function(modal, options, parameters) {
		this.modal = modal;
		this.element = $('#modal-form');
		this.options = options;
		this.parameters = parameters
		this.init();
	}

	DjangoModalAjaxForm.prototype.init = function() {
		this.initAjaxForm();
	}

	DjangoModalAjaxForm.prototype.initOnHideAfterSubmit = function(){
		var self = this;
		this.modal.on('hidden.bs.modal', function(){
			self.modal.remove();
			if(self.parameters.on_hide_modal_after_submit){
					self.parameters.on_hide_modal_after_submit();
			}
		});
	}

	DjangoModalAjaxForm.prototype.toogleSubmitState = function(){
		var new_state = (this.element.hasClass('disabled'))?'reset':'loading';
		$('#modal-submit').button(new_state)
	}

	DjangoModalAjaxForm.prototype.initAjaxForm = function(){
		var self = this;
		var n = false;
		var handler_response = new DjangoModalHandlerResponse();
		this.element.ajaxForm({
			dataType: 'json',
			data: self.options,
			success: function(response, statusText, xhr){
				var content = handler_response.handle(response);
				$('#modal-get-content').html(content);
				if(self.parameters.on_done){
					self.parameters.on_done(response, statusText, xhr);
				}
				self.initOnHideAfterSubmit();
				var newD = new DjangoModalAjaxForm(self.modal, self.options, self.parameters);
            },
			error: function(xhr, textStatus, errorThrown){
				self.initOnHideAfterSubmit();
				if(self.parameters.on_done){
					self.parameters.on_done(xhr, textStatus, errorThrown);
				}
                handler_response.handleErrorResponse(xhr, textStatus, errorThrown)
			},
			beforeSubmit: function(arr, form, df){
				self.toogleSubmitState();
				if(self.parameters.on_submit){
					self.parameters.on_submit();
				}
			},
		});

	}

	/*
		DjangoModal UtilRunner declaration
	*/

	var DjangoModalUtilRunner = function(modal, options, parameters){
		this.element = $('.util_runner');
		this.modal = modal;
		this.options = options;
		this.parameters = parameters;
		this.init();
	}

	DjangoModalUtilRunner.prototype.init = function(){
		var self = this;
		this.element.on('click', function(){
			self.toogleState();
			self.sendRequest();
			return false;
		});
	}
	DjangoModalUtilRunner.prototype.initOnHideAfterSubmit = function(){
		var self = this;
		this.modal.on('hidden.bs.modal', function(){
			self.modal.remove();
			if(self.parameters.on_hide_modal_after_submit){
					self.parameters.on_hide_modal_after_submit();
			}
		});
	}
	DjangoModalUtilRunner.prototype.sendRequest = function(){
		var handler_response = new DjangoModalHandlerResponse();
		var self = this;
		if(self.parameters.on_submit){
			self.parameters.on_submit();
		}
		var req = $.ajax({
			type: 'GET',
			dataType: 'json',
			url: self.element.attr('href'),
			data: self.options,
			success: function(response, statusText, xhr) {
                var content = handler_response.handle(response);
                $("#modal-get-content").html(content);
                //Re init the util runner because it's a new button
                var runner = new DjangoModalUtilRunner(self.modal, self.options, self.parameters);
			}, error: function(xhr, textStatus, errorThrown){
                handler_response.handleErrorResponse(xhr, textStatus, errorThrown)
			}, complete: function(){
				self.initOnHideAfterSubmit();
                if(self.parameters.on_done){
                    self.parameters.on_done();
                }
			},

		});
	}

	DjangoModalUtilRunner.prototype.toogleState = function(){
			var new_state = (this.element.hasClass('disabled'))?'reset':'loading';
			this.element.button(new_state)
	}

	/* 
		DjangoModal declaration
	*/
	var DjangoModal = function(modal_content, options, parameters){
		this.element = this.buildElement(modal_content);
		this.parameters = parameters;
		this.options = options;
		this.isUtil = this.hasUtilRunner();
		this.isForm = this.hasFormElement();
		this.init();
	}

	DjangoModal.prototype.init = function(){
		this.initOnHide();
		this.initOnShow();
		this.element.modal('show');
		if(this.isUtil){
			this.initUtilRunner();
		}
		if(this.isForm){
			this.initForm();
		}
	}

	DjangoModal.prototype.initUtilRunner = function(){
		var util_runner = new DjangoModalUtilRunner(this.element, this.options, this.parameters);
	}

	DjangoModal.prototype.initForm = function(){
		var form = new DjangoModalAjaxForm(this.element, this.options, this.parameters);
	}
	DjangoModal.prototype.initOnHide = function(){
		var self = this;
		this.element.on('hidden.bs.modal', function(){
			self.element.remove();
			if(self.parameters.on_hide_modal){
					self.parameters.on_hide_modal();
			}
		});
	}

	DjangoModal.prototype.initOnShow = function() {
		if(this.parameters.on_show_modal){
			this.element.on('show', this.parameters.on_show_modal());
		}
	}

	DjangoModal.prototype.buildElement = function(modal){
		$('body').prepend(modal);
		return $('#generic-modal');
	}

	DjangoModal.prototype.hasUtilRunner = function(){
		var res = false;
		if($('.util_runner').length){
			res = true;
		}
		return res;
	}

	DjangoModal.prototype.hasFormElement = function() {
		var res = false;
		if($('#modal-form').length){
			res = true;
		}		
		return res;
	}

	/*
		DjangoModalRunner declaration
	*/

	var DjangoModalRunner = function(element, modal_parameters){
		this.element = element;
		this.modal_parameters = modal_parameters;
		this.data_options = null;
		this.init();
	}
	DjangoModalRunner.prototype.init = function(){
		this.setOptions();
		//TODO init the loading modal
		this.sendRequest();
	}

	DjangoModalRunner.prototype.sendRequest = function(){
		var response_handler = new DjangoModalHandlerResponse();
		var self = this;
		var req = $.ajax({
			type: 'GET',
			dataType: 'json',
			url: this.element.attr('href'),
			data: self.data_options,
			success: function(response) {
				//TODO hide loading modal
				modal_content = response_handler.handle(response);
				var modal = new DjangoModal(modal_content, self.data_options, self.modal_parameters);
			}, error: function(jqXHR, textStatus, errorThrown){
				response_handler.handleErrorResponse(jqXHR, textStatus, errorThrown);
			}

		});

	}

	DjangoModalRunner.prototype.handleErrorResponse = function(text){

	}

	DjangoModalRunner.prototype.handleResponse = function(response){
		if(response.type === 'normal'){
			return this.handleNormalResponse(response);
		} else {
			return this.handleRedirectResponse(response);
		}
	}

	DjangoModalRunner.prototype.handleRedirectResponse = function(response){
		window.location = response.url;
	}

	DjangoModalRunner.prototype.handleNormalResponse = function(response){
		
	}

	DjangoModalRunner.prototype.setOptions = function() {
        var datas = {};
        var self = this;
        $.each(self.element.data(), function(key, value){
            if(self.element.attr('data-'+key) !== undefined){
                datas[key] = value;
            }
        });
		this.data_options = datas;
	}

	DjangoModalRunner.prototype.hasState = function(){
		var res = false;
		if(this.element.attr('data-loading-text')){
			res = true
		}
		return res;
	}

	$.fn.DjangoModalRunner=function(parameters)
    {
       var default_parameters = {
       		'events': ['click'],
            'on_event': null,
       		'on_show_modal': null,
       		'on_hide_modal': null,
       		'on_submit': null,
       		'on_hide_modal_after_submit': null,
       		'on_done': null,
       };

       var modal_parameters = $.fn.extend(default_parameters, parameters);
       return this.each(function()
       {
       		for(var i = 0; i < modal_parameters.events.length; i++){
       			$(this).on(modal_parameters.events[i], function(event){
                    if(modal_parameters.on_event){
                        modal_parameters.on_event(event);
                    }
           			new DjangoModalRunner($(this), modal_parameters);
           			return false;
       		    });
       		}
       });
    };

})( jQuery, window, document );
