  // Question is: do we take into account the no-synapses case? If
  // yes, we are forced to use MPI_Alltoallv; if no, we can use
  // MPI_Alltoall but some (how many?) packets are flying empty...

  // BEWARE: even if we are doing just one MPI_Alltoallv, there is
  // no way the receiver can expect a zero length packet, unless a
  // process already knows no synapses reach it, meaning
  // axSpikeForwardCount[h] == 0...

  // here I'm trying the Alltoall-with-empty-packets case

  // the fake_spike is signed with the max float value
  axonalSpikeDataOnlyClass fake_spike = {FLT_MAX, 0};
  for(h=0; h<lnp_par.globH; h++) {

    int lastRowInPacket = (h+1)*pktLength-1;

    // is 'pForwardAxonalSpikes[.].count' the correct count of how
    // many spikes I send to 'h'?
    fake_spike.spikes_in_packet = pForwardAxonalSpikes[h].count;
    axSpike_forwardBuffer[lastRowInPacket] = fake_spike;

    // spill 'pktLength-1' spikes into the packet, starting from the
    // tail and going backwards, then add a 'fake' last one with the
    // total number in the .spikes_in_packet field of the
    // 'pktLength'-th spike
    int count = (pForwardAxonalSpikes[h].count < pktLength
		 ? pForwardAxonalSpikes[h].count
		 : pktLength-1);

    // BEWARE: in non-full packets, we are leaving the unused
    // memory dirty

    // start of tail to pull spikes from (can be 0)
    int startTail = pForwardAxonalSpikes[h].count - count;

    // copy 'pktLength-1' spikes from the tail into
    // axSpike_forwardBuffer; they are contiguous both in the use
    // and the takeoff buffer, so memcpy should be equivalent but
    // a little faster than a 'for' loop
    memcpy(&axSpike_forwardBuffer[h*pktLength],
	   &pForwardAxonalSpikes[h].list[startTail],
	   count*sizeof(axonalSpikeDataOnlyClass));

    // we spilled some spikes, remember how many are left
    axSpike_forwardCount[h] = startTail;

    // if I want padding, here is the place to do it...

  }; // I've scanned each 'h' preparing fixed-length packets

  /* char filename3[20]; */
  /* sprintf(filename3, "dump_%u.pre_cpy_pre_ATA", lnp_par.loc_h); */

  /* save_and_tell(filename3, */
  /* 		(void*)axSpike_forwardBuffer, (void*)axSpike_backwardBuffer, */
  /* 		FC, BC, */
  /* 		FO, BO, */
  /* 		lnp_par.loc_h, */
  /* 		lnp_par.globH, */
  /* 		sizeof(typeof(axSpike_backwardBuffer[0]))); */

  MPI_Alltoall((void *) axSpike_forwardBuffer, pktLength, MPI_spike,
	       (void *) axSpike_backwardBuffer, pktLength, MPI_spike,
	       MPI_COMM_WORLD);


  /* char filename4[20]; */
  /* sprintf(filename4, "dump_%u.pre_cpy_post_ATA", lnp_par.loc_h); */

  /* save_and_tell(filename4, */
  /* 		(void*)axSpike_forwardBuffer, (void*)axSpike_backwardBuffer, */
  /* 		FC, BC, */
  /* 		FO, BO, */
  /* 		lnp_par.loc_h, */
  /* 		lnp_par.globH, */
  /* 		sizeof(typeof(axSpike_backwardBuffer[0]))); */

  /* // update the spikes count, to be able to dump... */
  /* for(h=0; h<lnp_par.globH; h++) */
  /*   FC[h]=pForwardAxonalSpikes[h].count; */

  /* // save the total spikes count, to be able to dump... */
  /* for(h=0; h<lnp_par.globH; h++) { */
  /*   int lastRowInPacket = (h+1)*pktLength-1; */
  /*   BC[h]=axSpike_backwardBuffer[lastRowInPacket].spikes_in_packet; */
  /* }; */


  // here I must decode each received packet's last_spike's to prepare
  // the following (eventually empty) MPI_Alltoallv

  for(h=0; h<lnp_par.globH; h++) {
    int lastRowInPacket = (h+1)*pktLength-1;

    // pull the real number of spikes from the last position in
    // the packet and save it in its final destination
    pBackwardAxonalSpikes[h].count =
      axSpike_backwardBuffer[lastRowInPacket].spikes_in_packet;

    // copy back as many spikes there are in the packet
    int count = (pBackwardAxonalSpikes[h].count < pktLength
		 ? pBackwardAxonalSpikes[h].count
		 : pktLength-1);

    // start of tail to attach spikes to (can be 0)
    int startTail = pBackwardAxonalSpikes[h].count - count;

    // copy 'pktLength-1' spikes from axSpike_backwardBuffer into
    // the tail; spikes are contiguous both in the landing buffer
    // and in their use buffer, so memcpy should be equivalent but
    // a little faster than a 'for' loop
    memcpy(&pBackwardAxonalSpikes[h].list[startTail],
	   &axSpike_backwardBuffer[h*pktLength],
	   count*sizeof(axonalSpikeDataOnlyClass));

    // we filled the tail with spikes, how many yet to fill
    axSpike_backwardCount[h] = startTail;

    // 'pBackwardAxonalSpikes[.].count' is set to the real size while
    // '.expectedCount' is unused, they are checked for equality later
    // on, so the former is dumped into the latter
    pBackwardAxonalSpikes[h].expectedCount = pBackwardAxonalSpikes[h].count;

  }; // I've scanned each 'h'

  /* char filename5[20]; */
  /* sprintf(filename5, "dump_%u.post_cpy", lnp_par.loc_h); */

  /* save_and_tell(filename5, */
  /* 		(void*)pForwardAxonalSpikes[0].list, (void*)pBackwardAxonalSpikes[0].list, */
  /* 		FC, BC, */
  /* 		axSpike_forwardOffset, axSpike_backwardOffset, */
  /* 		lnp_par.loc_h, */
  /* 		lnp_par.globH, */
  /* 		sizeof(typeof(pForwardAxonalSpikes[0].list[0]))); */
