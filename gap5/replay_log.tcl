# This script reads a gap5 .log file (written if DO_LOGGING is enabled) and
# attempts to produce a script that can be pasted into the gap5 console
# to reproduce the events recorded in the log.  It's not very elegant and
# is incomplete (e.g. disassemble reads is not supported) so the results may
# diverge from reality after a while.  It was useful for tracking down bugs
# with pair_rec updates and shuffle_pads, though.

# To use, it needs a copy of the database as it was before the events being
# replayed happened, and a .log file.  The script converts the .log file
# into a fragment of tcl script:

# tclsh replay_log.tcl < db.0.log > cmds

# Then open the database in gap5 and bring up the console window.  
# To ensure the right libraries are present, type 'edit_contig'.
# Then cut and paste the contents of the generated file into the console.

# N.B. As the replay may diverge from the original edits, it's possible
# for things to happen that wouldn't be allowed normally.  For example,
# if a break happens in a location slightly different from the original
# it's possible to clip sequences from the wrong contig editor later on,
# resulting in inconsistencies or worse when the contigs are saved.
# It's important to check for a faulty replay when finding the cause of
# database inconsistencies and crashes, to avoid wasting time on false
# leads.

# Map original contig numbers to new ones.  Needed because break_contig
# may return a different number to the one in the original log file.
proc map_ctg { ctg ccmap } {
    upvar $ccmap cmap
    if { [ info exists cmap($ctg) ] } {
	return $cmap($ctg)
    }
    return $ctg
}

set last_pid_host {}
set cmds {}
set main_io {}
array set saved {}
array set join_eds {}
array set contig_editors {}
array set editor_contigs {}

# Parse the log file.
while {[gets stdin line] >= 0} {

    # Ensure we only replay data from a single gap5 run by checking
    # that the host name and pid don't change.
    regexp {\[\d+@\w+\]} $line pid_host
    if { $last_pid_host eq "" } { set last_pid_host $pid_host }
    if { $last_pid_host ne $pid_host } { break }

    if { [regexp {log_call contig_editor (\.e\d+) -io io=(0x[0-9a-f]+) -contig (\d+)} $line cmd ed_num io ctg_num ] } {

	# Parse contig_editor command.  Keep track of which editors
	# are open on each contig, and which contig is in each editor.

	if { ! [regexp { \-pos (-?\d+)} $line dummy pos] } { set pos "" }
	if {[regexp { \-contig2 (\d+)} $line dummy ctg2]} {
	    # Join editor version
	    if { ! [regexp { \-pos2 (-?\d+)} $line dummy pos2] } { set pos2 "" }
	    lappend cmds [list contig_editor $ed_num $ctg_num $ctg2 $pos $pos2]
	    array set join_eds [list "$ctg_num/$ctg2" $ed_num]
	    lappend contig_editors($ctg_num) "$ed_num.ed1.pane.seq.sheet"
	    lappend contig_editors($ctg2) "$ed_num.ed0.pane.seq.sheet"
	    lappend editor_contigs($ed_num) $ctg_num $ctg2
	} else {
	    # Ordinary contig editor
	    lappend cmds [list contig_editor $ed_num $ctg_num $pos]
	    lappend contig_editors($ctg_num) "$ed_num.ed1.pane.seq.sheet"
	    lappend editor_contigs($ed_num) $ctg_num
	}
	puts stderr ":: $ed_num $editor_contigs($ed_num)"
	if { $main_io eq "" } {
	    set main_io $io
	}
	array set saved [list $ed_num -1]

    } elseif { [regexp { log_str editor_exit (\.e\d+) } $line cmd ed_num ] } {

	# Parse editor_exit

	lappend cmds [list editor_exit $ed_num ]
	puts stderr "exit $ed_num"
	if { [info exists editor_contigs($ed_num)] } {
	    foreach ctg $editor_contigs($ed_num) {
		set new_list {}
		foreach ed $contig_editors($ctg) {
		    # puts stderr "xx $ctg $ed"
		    if { [regexp {^(\.e\d+)\.} $ed dummy e] } {
			# puts stderr "xx $ctg $ed $e $ed_num"
			if { $e ne $ed_num } {
			    lappend new_list $ed
			}
		    }
		}
		puts stderr "yy $ctg $new_list"
		set contig_editors($ctg) $new_list
		set editor_contigs($ed_num) {}
	    }
	}

    } elseif { [regexp { log_call (\.e\d+)\.ed(\d)\.pane\.seq\.sheet save} $line cmd ed_num ed] } {

	# Parse contig editor save button

	set s [array get saved $ed_num]
	if { $s eq "" } {
	    puts stderr "save on unknown editor $ed_num"
	    exit 1
	}
	array set saved [list $ed_num [llength $cmds]]

	# If the previous command was editor_exit, user closed the window
	# then said 'yes' to save.  In this case we need to reverse the
	# order of events to we save before exit.

	set prev [lindex $cmds end]
	if {[lindex $prev 0] eq "editor_exit" && [lindex $prev 1] eq $ed_num} {
	    lset cmds end [list save $ed_num $ed]
	    lappend cmds $prev
	} else {
	    lappend cmds [list save $ed_num $ed]
	}

    } elseif { [regexp { log_call (\.e\d+)\.ed(\d)\.pane\.seq\.sheet join} $line cmd ed_num ed] } {

	# Parse join editor join button

	array set saved [list $ed_num [llength $cmds]]
	lappend cmds [list join $ed_num $ed]

	foreach { ctg1 ctg2 } $editor_contigs($ed_num) break
	foreach ed $contig_editors($ctg2) {
	    if { [regexp {^(\.e\d+)\.} $ed dummy e] } {
		lappend contig_editors($ctg1) "$ed"
		set editor_contigs($e) $ctg1
	    }
	}
	set contig_editors($ctg2) {}

    } elseif { [regexp { edJoinAlign fixed_left=(\d) fixed_right=(\d) =(\d+)@(-?\d+) =(\d+)@(-?\d+)} $line cmd fleft fright ctg1 pos1 ctg2 pos2] } {

	# Parse join editor align button

	set ed [array get join_eds "$ctg2/$ctg1"]
	if {$ed ne ""} {
	    lappend cmds [list edJoinAlign [lindex $ed 1] $fleft $fright $ctg1 $pos1 $ctg2 $pos2]
	} else {
	    puts stderr "No join editor for $ctg1 $ctg2"
	    exit 1
	}

    } elseif { [regexp { editor_delete_base (\.e\d+)\.ed(\d)\.pane\.seq\.sheet {(\d+) (\d+) (-?\d+)} (\d) (\d) (\d)} $line cmd ed_num ed type ctg pos end dir powerup] } {

	# Parse editor_delete_base

	lappend cmds [list editor_delete_base $ed_num $ed $type $ctg $pos $end $dir $powerup]

    } elseif { [regexp { editor_clip_seq (\.e\d+)\.ed(\d)\.pane\.seq\.sheet {(\d+) (\d+) (-?\d+)} ([lr])} $line cmd ed_num ed type ctg pos dirn] } {

	# Parse editor_clip_seq

	lappend cmds [list editor_clip_seq $ed_num $ed $type $ctg $pos $dirn]

    } elseif { [regexp { editor_clip_contig (\.e\d+)\.ed(\d)\.pane\.seq\.sheet {(\d+) (\d+) (-?\d+)} ([lr])} $line cmd ed_num ed type ctg pos dirn] } {

	# Parse editor_clip_contig
	# Only bother to record if it wasn't preceeded by editor_clip_seq
	# on the same contig.

	if { [llength $cmds] > 0 \
		 && !( [lindex $cmds end 0] eq "editor_clip_seq"
		      && [lindex $cmds end 1] eq "$ed_num"
		      && [lindex $cmds end 3] == $type
		      && [lindex $cmds end 4] == $ctg) } {
	    lappend cmds [list editor_clip_contig $ed_num $ed $type $ctg $pos $dirn]
	}

    } elseif { [regexp { break_contig -io io=0x[0-9a-f]+ -contig (\d+) -pos (-?\d+) -break_holes (\d)} $line cmd ctg pos holes] } {

	# Parse break_contig

	lappend cmds [list break_contig 0 $ctg $pos $holes]

    } elseif { [regexp {shuffle_pads -io io=(0x[0-9a-f]+) -contigs {((?:{=\d+ -?\d+ -?\d+})+)} -flush (\d+) -band (\d+) -soft_clips (\d+) -max_pass (\d+)} $line cmd io contigs do_flush band soft_clips max_pass] } {

	# Parse shuffle_pads

	if { $io eq $main_io } {
	    # Verision invoked from the main menu
	    lappend cmds [list shuffle_pads main $contigs $do_flush $band $soft_clips $max_pass]
	} else {
	    # Version invoked from a contig editor.  We need to work out
	    # which one it was, as this isn't recorded in the log file.
	    puts stderr "$contigs"
	    foreach { item } $contigs {
		puts stderr "$item"
		foreach { eqctg from to } $item { 
		    set ctg [string range $eqctg 1 end]
		    puts stderr "$eqctg $ctg $contig_editors($ctg)"
		    if { [llength $contig_editors($ctg)] > 0 } {
			set ed [lindex $contig_editors($ctg) 0]
			lappend cmds [list shuffle_pads $ed [list [list $ctg $from $to]] $do_flush $band $soft_clips $max_pass]
		    } else {
			lappend cmds [list shuffle_pads {unknown} [list [list $ctg $from $to]] $do_flush $band $soft_clips $max_pass]
		    }
		}
	    }
	}

    } elseif {[regexp { io_undo_exec (\.e\d+)\.ed(\d)\.pane\.seq\.sheet (\d+)} $line cmd ed_num ed ctg]} {

	# Parse io_undo_exec

	lappend cmds [list io_undo $ed_num $ed $ctg]
    }
}


# Generate the new script fragment.

set p 0
array set used {}
array set cmap {}
set new_ctg 0
for { set ci 0 } { $ci < [llength $cmds] } { incr ci } {
    set cmd [lindex $cmds $ci]
    set c [lindex $cmd 0]
    set ed_num [lindex $cmd 1]
    set s [array get saved $ed_num]
    if { $s ne "" } {
	set sp [lindex $s 1]
    } else {
	set sp -1
    }
    puts "# $cmd"
    if { 1 || $p <= $sp \
	     || ($c eq "editor_exit" && [array get used $ed_num] ne "") \
	     || $c eq "break_contig" } {
	if { $c eq "contig_editor" } {
	    array set used [list $ed_num 1]
	}

	if { $c eq "contig_editor" } {

	    # contig_editor command

	    if {[llength $cmd] == 4} {
		# Ordinary contig editor
		set ctg1 [map_ctg [lindex $cmd 2] cmap]
		if {[lindex $cmd 3] ne ""} {
		    set pos " -pos [lindex $cmd 3]"
		} else {
		    set pos ""
		}
		puts "contig_editor $ed_num -io \$io -contig $ctg1$pos"
	    } else {
		# Join editor
		set ctg1 [map_ctg [lindex $cmd 2] cmap]
		set ctg2 [map_ctg [lindex $cmd 3] cmap]
		if {[lindex $cmd 4] ne ""} {
		    set pos " -pos [lindex $cmd 4]"
		} else {
		    set pos ""
		}
		if {[lindex $cmd 5] ne ""} {
		    set pos2 " -pos2 [lindex $cmd 5]"
		} else {
		    set pos2 ""
		}
		puts "contig_editor $ed_num -io \$io -contig $ctg1$pos -contig2 $ctg2$pos2"
	    }
	    # Call update so the editor appears
	    puts "update idletasks"

	} elseif { $c eq "editor_exit" } {

	    # editor_exit command

	    puts "editor_exit $ed_num 0 1"
	    puts "update idletasks"

	} elseif { $c eq "save" } {

	    # contig editor save

	    set ed [lindex $cmd 2]
	    puts "log_call $ed_num.ed$ed.pane.seq.sheet save"

	} elseif { $c eq "join" } {

	    # contig editor join

	    set ed [lindex $cmd 2]
	    puts "$ed_num.ed$ed.pane.seq.sheet join"

	} elseif { $c eq "editor_delete_base" } {

	    # editor_delete_base
	    # Need to set cursor to the right position first for this to
	    # work properly.

	    foreach { ed type ctg pos end dir powerup } [lrange $cmd 2 end] {
		set mctg [map_ctg $ctg cmap]
		puts "$ed_num.ed$ed.pane.seq.sheet set_cursor 17 $mctg $pos"
		puts "editor_delete_base $ed_num.ed$ed.pane.seq.sheet \[list $type $mctg $pos\] $end $dir $powerup"
		puts "update idletasks"
	    }

	} elseif { $c eq "editor_clip_seq" } {

	    # editor_clip_seq

	    foreach { ed type ctg pos dirn } [lrange $cmd 2 end] {
		set mctg [map_ctg $ctg cmap]
		puts "editor_clip_seq $ed_num.ed$ed.pane.seq.sheet \[list $type $mctg $pos\] $dirn"
	    }

	} elseif { $c eq "editor_clip_contig" } {

	    # editor_clip_contig

	    foreach { ed type ctg pos dirn } [lrange $cmd 2 end] {
		set mctg [map_ctg $ctg cmap]
		puts "editor_clip_contig $ed_num.ed$ed.pane.seq.sheet \[list $type $mctg $pos\] $dirn"
	    }

	} elseif { $c eq "edJoinAlign" } {

	    # contig editor align button

	    foreach { fleft fright ctg1 pos1 ctg2 pos2 } [lrange $cmd 2 end] {
		puts "set $ed_num\(Lock) 0"
		puts "editor_lock $ed_num"
		puts "$ed_num.ed1.pane.seq.sheet xview $pos1"
		puts "$ed_num.ed0.pane.seq.sheet xview $pos2"
		puts "set $ed_num\(Lock) 1"
		puts "editor_lock $ed_num"
		puts "$ed_num.ed0.pane.seq.sheet join_align $fleft $fright"
	    }

	} elseif { $c eq "break_contig" } {

	    # break_contig
	    # We have to capture the record number for the new contig,
	    # and add it to cmap so we can translate the number that was
	    # used in the log file.  We get the old contig number by
	    # looking for the second contig editor to come up after
	    # break contig.

	    foreach { ctg pos holes } [lrange $cmd 2 end] {
		set mctg [map_ctg $ctg cmap]
		incr new_ctg
		puts "set new_ctg_$new_ctg \[break_contig -io \$io -contig $mctg -pos $pos -break_holes $holes\]"
	    }
	    if { [lindex $cmds [expr { $ci + 2 } ] 0 ] eq "contig_editor" } {
		set ctg [lindex $cmds [expr { $ci + 2 } ] 2]
		array set cmap [list $ctg "\$new_ctg_$new_ctg"]
	    } else {
		puts stderr "Couldn't get old contig name for break_contig"
		exit 1
	    }

	} elseif { $c eq "shuffle_pads" } {

	    # Shuffle pads

	    foreach { ed contigs do_flush band soft_clips max_pass } [lrange $cmd 1 end] {
		if { $ed eq "unknown" } {
		    # puts -nonewline "# "
		}
		puts -nonewline "log_call shuffle_pads -io "
		if { $ed eq "main" || $ed eq "unknown" } {
		    # Shuffle pads from main menu
		    puts -nonewline "\$io -contigs \[list"
		} else {
		    # Shuffle pads from contig editor
		    puts -nonewline "\[$ed io\] -contigs \[list"
		}
		# Translate the contig numbers if necessary
		foreach { contig } $contigs {
		    foreach { ctg from to } $contig break
		    set mctg [map_ctg $ctg cmap]
		    puts -nonewline " \[list \"=${mctg}\" $from $to\]"
		}
		puts "\] -flush $do_flush -band $band -soft_clips $soft_clips -max_pass $max_pass"
	    }

	} elseif { $c eq "io_undo" } {

	    # undo pressed.

	    foreach { ed ctg } [lrange $cmd 2 end] break
	    set mctg [map_ctg $ctg cmap]
	    puts "io_undo $ed_num.ed$ed.pane.seq.sheet $mctg"
	}
    }
    incr p
}
