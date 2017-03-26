import execjs

if1 = open('jsme.nocache.js','r').readlines()

execjs.eval(if1)

ctx = execjs.compile("""

    function jsmeOnLoad() {
        jsmeApplet = new JSApplet.JSME("jsme_container", "380px", "340px");
   }
   
   """)
   

