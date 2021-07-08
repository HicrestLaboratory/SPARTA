# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 11:12:55 2021

@author: Paolo
"""



def check_fixed_condition(val, condition):
    if type(condition) is list:
        return (val >= condition[0] and val <= condition[1]);
    else:
        return val == condition;
    
def check_good_sample(i):
    return (experiments["total_nonzeros"][i] != 0);


def plot_x_vs_y(x_field_function, y_field_function, fixed, label = None, draw_std_error = True):
    y_lists = {};
    y_values = [];
    y_errors =[];
    for i in range(n_exp):

        #check condition for adding experiment
        skip = False; 
        for fixed_field, condition in fixed.items():
            if not check_fixed_condition(experiments[fixed_field][i], condition):
                skip = True;
        if not check_good_sample(i):
            skip = True;
        if skip:
            continue;   
        
        x = x_field_function(i);
        if x not in y_lists:
            y_lists[x] = [];
        
        y = y_field_function(i);
        
        #append to data list
        y_lists[x].append(y);
        
    for x in y_lists.keys():
        y_values.append(np.mean(y_lists[x]));
        y_errors.append(np.std(y_lists[x]));

    x_values = list(y_lists.keys());
    plt.scatter(x_values,y_values, label = label);
    
    if draw_std_error:
        plt.fill_between(x_values, np.array(y_values) - np.array(y_errors),
                             np.array(y_values) + np.array(y_errors), alpha=0.1,
                             color="r");

def plot_x_vs_y_all_z(x_field_function, y_field_function, fixed, varying, draw_std_error = True):
    
    varying_values = find_unique(varying);
    for val in varying_values:
        new_fixed = dict(fixed);
        new_fixed[varying] = val;
        plot_x_vs_y(x_field_function, y_field_function, new_fixed, label = val);
    
    return 0

def make_title(fixed_values_dict):
    title = "";
    for key, value in fixed_values_dict.items():
        title += key + ":" + str(value) + "\n"
    return title
