// $Id$

/**
 * @file
 * Drives our implementation 
 * of the checkbox tree using dynatree 
*/


/**
 * Prepare interactive tree for gmod_dbsf_phylogeny_tree
 */
 $(function(){
      $("#checktree").dynatree({
        //Tree parameters
        persist: false,
        checkbox: true,
        selectMode: 3,
        activeVisible: true,
  
         onDblClick: function(dtnode, event) {
        dtnode.toggleSelect();
        },
        onKeydown: function(dtnode, event) {
         if( event.which == 32 ) {
          dtnode.toggleSelect();
          return false;
         }
        },

        //Un/check real checkboxes recursively after selection
        onSelect: function(select, dtnode) {
          dtnode.visit(function(dtnode){
            $("#edit-"+dtnode.data.key).attr("checked",select);
            },null,true);
        },
        //Hack to prevent appearing of checkbox when node is expanded/collapsed
        onExpand: function(select, dtnode) {
          $("#edit-"+dtnode.data.key).attr("checked",dtnode.isSelected()).addClass("hidden-checkbox");
        }
      });
      //Update real checkboxes according to selections
      $.map($("#checktree").dynatree("getTree").getSelectedNodes(),
       function(dtnode){
          $("#edit-"+dtnode.data.key).attr("checked",true);
          dtnode.activate();
        });
      });

