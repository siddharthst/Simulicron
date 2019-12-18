# Refer to:
# https://stackoverflow.com/questions/1389180/
###  StdLib  ###
from functools import wraps
import inspect

###########################################################################################################################
#//////|   Decorator   |//////////////////////////////////////////////////////////////////////////////////////////////////#
###########################################################################################################################

def auto_assign_arguments(function):

  @wraps(function)
  def wrapped(self, *args, **kwargs):
    _assign_args(self, list(args), kwargs, function)
    function(self, *args, **kwargs)

  return wrapped


###########################################################################################################################
#//////|   Utils   |//////////////////////////////////////////////////////////////////////////////////////////////////////#
###########################################################################################################################

def _assign_args(instance, args, kwargs, function):

  def set_attribute(instance, parameter, default_arg):
    if not(parameter.startswith("_")):
      setattr(instance, parameter, default_arg)

  def assign_keyword_defaults(parameters, defaults):
    for parameter, default_arg in zip(reversed(parameters), reversed(defaults)):
      set_attribute(instance, parameter, default_arg)

  def assign_positional_args(parameters, args):
    for parameter, arg in zip(parameters, args.copy()):
      set_attribute(instance, parameter, arg)
      args.remove(arg)

  def assign_keyword_args(kwargs):
    for parameter, arg in kwargs.items():
      set_attribute(instance, parameter, arg)
  def assign_keyword_only_defaults(defaults):
    return assign_keyword_args(defaults)

  def assign_variable_args(parameter, args):
    set_attribute(instance, parameter, args)

  POSITIONAL_PARAMS, VARIABLE_PARAM, _, KEYWORD_DEFAULTS, _, KEYWORD_ONLY_DEFAULTS, _ = inspect.getfullargspec(function)
  POSITIONAL_PARAMS = POSITIONAL_PARAMS[1:] # remove 'self'

  if(KEYWORD_DEFAULTS     ): assign_keyword_defaults     (parameters=POSITIONAL_PARAMS,  defaults=KEYWORD_DEFAULTS)
  if(KEYWORD_ONLY_DEFAULTS): assign_keyword_only_defaults(defaults=KEYWORD_ONLY_DEFAULTS                          )
  if(args                 ): assign_positional_args      (parameters=POSITIONAL_PARAMS,  args=args                )
  if(kwargs               ): assign_keyword_args         (kwargs=kwargs                                           )
  if(VARIABLE_PARAM       ): assign_variable_args        (parameter=VARIABLE_PARAM,      args=args                )
###########################################################################################################################