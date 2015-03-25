CodeReview
##########

This is the code I used as part of a monthly code review session in the [GJB group](http://www.compbio.dundee.ac.uk/)


### Description
I took some bits and pieces from my framework codebase to make the program below. 
This particular program gets variation data (either variants or mutations) for a given UniProt identifier (or multiple).
Regarding the code itself, some comments or doc strings might be wrong because I took some bits out for making it shorter. 
The program contains a Sequence Class which holds sequence information (focusing on variation data in this case) and
a command line parser. As you can see immediately on the top, main.py loads a number of parsers, fetchers and 
utils which are used in the class methods. Those are not for review but check them as well if you are curious…

### Usage

First you need to checkout the code to your cwd. 

A simple example of use if you want to see the all the steps and JSON output, would be to run `python main.py -v -i P00439` 


otherwise simply run (it takes a bit to complete if the number of variants is big) `python main.py -i P00439 > output.json`

### Dependencies
python 2.7+

requests - `pip install requests`

If you want to try out the ipython notebook provided ('Code Review - Example.ipynb') 
you just need to install `pip install ipython[notebook]`

To run the ipython notebook just cd to your ~CodeReview/ and run `ipython notebook`


### Some of the comments I got


""
very nice, cleanly laid out and commented
could maybe decompose a few methods further to keep compact, e.g. init to loadFasta(), loadText(), or the 
JSON parse loop in fetch_variants_from_ensembl_rest
utils.py:
request_info_url could maybe also check status code 301/302 (redirects), for industrial strength
fetchers.py:
is line.rstrip(“\r\n”) platform portable (do you care)?
test at 'filtering out synonymous variants' has a very complex condition - hard to maintain/debug - could
 expand with temporary variables for readability
these are all picky comments, nothing major at all!
"""

"""
You are misusing Try/Except. If you are not fishing for an specific error, use if/else:
For instance, it would be more readable

if len(identifier) != 6:

  message = "Error: UniProt identifiers are usually 6-character long..."
  if self.verbose:
      log(message)
  raise ValueError(message)
* I can't understand the usage of self.uniprot_txt. 
* JSON scheme support integers, it would be better to output sequence positions as integers.
"""

"""
There might be something up with your summary section of your JSON?

  "SUMMARY": {
      "UNIPROT_NAME": "MSTAVLENPGLGRKLSDFGQETSYIEDNCNQNGAISLIFSLKEEVGALAKVLRLFEENDVNLTHIESRPSRLKKDEYEFFTHLDKR
      SLPALTNIIKILRHDIGATVHELSRDKKKDTVPWFPRTIQELDRFANQILSYGAELDADHPGFKDPVYRARRKQFADIAYNYRHGQPIPRVEYMEEEKKTWGT
      VFKTLKSLYKTHACYEYNHIFPLLEKYCGFHEDNIPQLEDVSQFLQTCTGFRLRPVAGLLSSRDFLGGLAFRVFHCTQYIRHGSKPMYTPEPDICHELLGHVP
      LFSDRSFAQFSQEIGLASLGAPDEYIEKLATIYWFTVEFGLCKQGDSIKAYGAGLLSSFGELQYCLSEKPKLLPLELEKTAIQNYTVTEFQPLYYVAESFND
      AKEKVRNFAATIPRPFSVRYDPYTQRIEVLDNTQQLKILADSINSEIGILCSALQKIK" .......
      
I'm fairly certain that's not the Uniprot name? Also: Remind me to ask what the flash method does.
"""

"""
If you require version dependant packages its a good idea to catch the errors related to this and provide a 
sensible fail that suggests that the version might be a problem, or better yet check the version of python 
in the script at the start and fail gracefully if it is insufficient.
Only really major problem is the use of catch-all try/excepts - these should not be used (except in very 
extreme cases) - try/except is for catching specific erros, not replacing if/else. Thigo has the right of it.
Also, main_handler is kind of redundant. You can call class methods from the constructor (although I don't 
know if this is considered bad style) so if these are going to be called all the time, every time, I'd do it there.
No check on the input before it gets used as a string and then the value is checked to make sure its 6 
chars long. What about if the input given doesn't have a string method (yes, it should be hey), and what 
happens if you give it, say, int(123456)? Answer: “Loading UNIPROT ID 123456…, Error: UniProt identifiers 
are usually 6-character long…” which is not a useful error in this case since the problem is that it isn't a 
string rather than that it isn't 6 chars long. What about or [“P”,”0”,”0”,”4”,”3”,”9”]? “Loading 
UNIPROT ID ['P', '0', '0', '4', '3', '9']… 404 http://www.uniprot.org/uniprot/['P', '0', '0', '4', '3', '9'].fasta 
Warning: ['P', '0', '0', '4', '3', '9'].fasta not available for download.” Would be good practice to check 
that the input can be cast as a sensible-looking string before using it.
"""
