---
title : Creating New Test Cases
layout : page_getting-started
permalink: /fehmpytest-info/newtest
hero_height: is-hidden
---

# Creating New Test Cases

### Topics

* [Installation](/FEHM/getting-started/fehmpytest-info/install)
* [Testing FEHM](/FEHM/getting-started/fehmpytest-info/testing)
* Creating New Test-Cases
* [Test-Case Description](/FEHM/getting-started/fehmpytest-info/testdesc)

---


A developer can add new test-cases to the suite. There are two steps 
to adding new test-cases: 

    1. Create the test-case folder.
    2. Add the test method.

## Create the Test-Case Folder


This is the folder structure of a test-case::

    [test-case]
        |
        |_[input]
        |   |
        |   |_[control]
        |  
        |_[compare]

All input files needed to run the FEHM functionality of the test-case go inside
the input folder. All control files go inside the control folder. All compare 
files (contour, history, and output) that are known to be correct go inside the
compare folder.

To set up a new test-case folder:

1. Go into the folder *fehmpytests*.
2. Create a folder ```<test-case>``` where ```<test-case>``` is the name of the new 
   case.
3. Inside ```<test-case>```, create two folders called *input* and *compare*.
4. Inside the *input* folder, create a folder called *control*.
5. In the *control* folder, place all control files.
6. If there is only one control file, rename it to **fehmn.files**.
7. If there are more than one control file, rename each file to 
   ```<subcase>.files``` where ```<subcase>``` is the name of the subcase.
8. In the *input* folder, place all input files needed for the FEHM run.
9. In the *compare* folder, place all comparison files known to be correct.

     
## Add the Test Method


To add the test method:

1. Open **fehmpytests.py**.
2. Inside the class 'Tests', write a method ```test_<name>``` where ```<name>``` is 
   the name of the test-case. Here is an example for the *avdonin* test:

   ```
   
       class Tests(unittest.TestCase):
           
           ...
           
           #This is the new test method for avdonin.
           def test_avdonin(self):
               ...
               
           ...    
```   
3. Inside this test method, call 
       
   ``self.test_case('<test-case>')`` 
       
   where ``<test-case>`` is the name of the folder you created for the new 
   test-case. See ``fehmpytests.fehmTest.test_case`` below for details on the 
   general test case method. Here is an example for the test method for
   **avdonin**:
   ```
       class Tests(unittest.TestCase):
           
           ...
           
           def test_avdonin(self):
               #Add this line to call the general test case for avdonin.
               self.test_case('avdonin')
               
           ...
   ```         
4. Inside the class *Suite*, under the condition *all*, add the following line:
 ```
       suite.addTest(Tests('<method-name>'))
   ```    
   where ```<method-name>``` is the name of the test method you just defined. Here 
   is an example for adding the **avdonin** test to the test-suite:
   ```
       def suite(case, test_case):
           suite = unittest.TestSuite()
            
           if case == 'all':
               suite.addTest(Tests('test_saltvcon'))
               suite.addTest(Tests('test_dissolution'))
               suite.addTest(Tests('test_salt_perm_poro'))
               
               #This is how the avdonin test is added to the test-suite.
               suite.addTest(Tests('test_avdonin'))
               
               ...
 ```  
Running fehmpytests will now include the new test-case.     

## Customizing a Test-Case


By default, *test_case()* will check for a maximum difference less than 1.e-4
on all attributes of the FEHM simulation. Passing a dictionary into 
*test_case()* as the second argument allows a developer to specify how these 
tests are performed. The following keywords are recognized by test_case():

    + 'variables':    list   
    + 'nodes':        list   
    + 'components':   list    
    + 'maxerr':       float   
    + 'test_measure': string 
    
The following is an example for specifying the components, variables, and format
for the **saltvcon** test:
 ```   
    #Pass a dictionary into test_case() with keywords specifed.
    def test_saltvcon(self):
        arguments = {}
        arguments['components'] = ['water']
        arguments['variables']  = ['Kx']
        arguments['format'] = 'relative' 
          
        self.test_case('saltvcon', arguments)  
 ```       
## Documentation on test_case() Method

```
fehm.test_case(name, parameters={})
```

Performs a test on a FEHM simulation and raises an AssertError if it fails the test.<br><br>


##### Parameters:

**name (str)** - The name of the test-case folder.<br>
**parameters (dict)** - Attribute values that override default values.

* Key Value Choices
	- 'variables': list[str, str, ...]
	- 'times': list[float, float, ...]
	- 'nodes': list[int, int, ...]
	- 'components': list[str, str, ...]
	- 'maxerr': float
	- 'test_measure': str
		* 'max_difference'
		* 'rms_difference'
		* 'perc_difference'                                  
            
The folder 'name' in fehmpytests must exist with correct structure.If parameters are not passed into this method, all simulated attributes will be checked for a relative difference of less than 1.e-4.

